# ==========================================================================================================
# Class for funGp models
# ==========================================================================================================

# Constructor of the class
# ----------------------------------------------------------------------------------------------------------
#' @title Class: functional input Gaussian process model
#' @description To create a funGp object, use \link[funGp]{funGp} . See also this function for mode details.
#'
#' @slot call Object of class \code{"language"}. User call reminder.
#' @slot type Object of class \code{"character"}. Type of model based on inputs structure. To be chosen from {"scalar", "functional", "hybrid"}.
#' @slot ds Object of class \code{"numeric"}. Number of scalar inputs.
#' @slot df Object of class \code{"numeric"}. Number of functional inputs.
#' @slot f_dims Object of class \code{"numeric"}. Dimension of each functional input.
#' @slot sIn Object of class \code{"matrix"}. Scalar inputs.
#' @slot fIn Object of class \code{"list"}. Functional inputs. Each element of the list contains a functional input in the form of a matrix.
#' @slot sOut Object of class \code{"matrix"}. Scalar output.
#' @slot n.tot Object of class \code{"integer"}. Number of observed points (not necessarily all are used for prediction).
#' @slot n.tr Object of class \code{"integer"}. Number of training points.
#' @slot f_proj Object of class \code{"funGpProj"}. Data structures related to the projection of functional inputs.
#' @slot kern Object of class \code{"funGpKern"}. Data structures related to the kernel.
#' @slot preMats Object of class \code{"list"}. L and LInvY matrices pre-computed for prediction. L is a lower diagonal matrix such that
#' \eqn{L'L} equals the training cross covariance matrix \eqn{K.tt}. On the other hand, \eqn{LInvY = L^(-1) * sOut}.
#'
#' @rdname funGp-class
#' @include 2_funGpProj_Class.R
#' @include 2_funGpKern_Class.R
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setClass("funGp",
         representation(
           call = "language",          # user call reminder
           type = "character",         # Type of model. To be chosen from {"scalar", "functional", "hybrid"}.
           ds = "numeric",             # number of scalar inputs
           df = "numeric",             # number of functional inputs
           f_dims = "numeric",         # dimension of each functional input
           sIn = "matrix",             # scalar inputs
           fIn = "list",               # each element (n x fDims_i) contains a functional input
           sOut = "matrix",            # scalar output
           n.tot = "integer",          # number of observed points
           n.tr = "integer",           # number of training points
           f_proj = "funGpProj",         # structures related to the projection of functional inputs
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
#' reduction.
#' @param kerType an optional character specifying the covariance structure to be used. To be chosen between
#' "gauss", "matern5_2" and "matern3_2". Default is "matern5_2".
#' @param f_disType an optional character specifying the distance function to use for the functional inputs
#' within the covariance function. To be chosen between "scalar" and "functional". Default is "functional".
#' @param f_pdims an optional array with the projection dimension for each functional input.
#' @param f_family fill!!!!!!
#' @param var.hyp fill!!!!!!
#' @param ls_s.hyp fill!!!!!!
#' @param ls_f.hyp fill!!!!!!
#' @param n.starts Fill!!!!!!!!!!
#' @param n.presample Fill!!!!!!!!!!
#'
#' @importFrom methods new
#' @importFrom stats optim
#' @importFrom stats runif
#'
#' @examples
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), matrix(runif(n.tr*22), ncol = 22))
#'
#' # generating output data for training
#' sOut <- as.matrix(sapply(t(1:n.tr), function(i){
#'   x1 <- sIn[i,1]
#'   x2 <- sIn[i,2]
#'   f1 <- fIn[[1]][i,]
#'   f2 <- fIn[[2]][i,]
#'   t1 <- seq(0,1,length = length(f1))
#'   t2 <- seq(0,1,length = length(f2))
#'   as.numeric(x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2))
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
funGp <- function(sIn = NULL, fIn = NULL, sOut, kerType = "matern5_2",
                  f_disType = "functional", f_pdims = 3, f_family = "PCA",
                  var.hyp = NULL, ls_s.hyp = NULL, ls_f.hyp = NULL,
                  n.starts = 1, n.presample = 20) {
  # extend simplified user inputs to full versions
  if (!is.null(fIn)) {
    if (length(f_disType) == 1) f_disType <- rep(f_disType, length(fIn))
    if (length(f_pdims) == 1) f_pdims <- rep(f_pdims, length(fIn))
    if (length(f_family) == 1) f_family <- rep(f_family, length(fIn))
  }

  # check validity of user inputs
  checkVal_funGp(as.list(environment()))

  # create objects of class funGpKern and funGp
  kern <- new("funGpKern")
  model <- new("funGp")

  # extract generic information from user inputs
  sOut <- as.matrix(sOut)
  n.tr <- length(sOut)
  n.presample <- max(n.presample, n.starts)

  # 3 possible cases
  # Case 1: scalar and functional
  # Case 2: functional only
  # Case 3: scalar only
  if (all(!is.null(sIn), !is.null(fIn))) { # Hybrid-input case *******************************************
    # extract information from user inputs specific to the hybrid-input case
    sIn <- as.matrix(sIn)
    ds <- ncol(sIn)
    df <- length(fIn)
    f_dims <- sapply(fIn, ncol)

    # perform projection of functional inputs
    # the projection is such that F = X * B' + e, with
    #
    # n: input points
    # k: original dimension
    # p: projection dimension
    # F: original inputs ............. matrix of dimension nxk
    # B: basis functions ............. matrix of dimension pxk (one basis per column)
    # X: projection coefficients ..... matrix of dimension nxp
    bcj <- dimReduction(fIn, df, f_pdims, f_family)
    f_basis <- bcj$basis
    f_coefs <- bcj$coefs
    f_J <- bcj$J

    # compute scalar distance matrices
    sMs <- setScalDistance(sIn, sIn)

    # compute functional distance matrices
    fMs <- setFunDistance(f_coefs, f_coefs, f_J)

    # optimize hyperparameters if some is required
    if (all(!is.null(var.hyp), !is.null(ls_s.hyp), !is.null(ls_f.hyp))) {
      varHyp <- var.hyp
      lsHyps <- c(ls_s.hyp, ls_f.hyp)
    } else {
      hypers <- setHypers_SF(sMs, fMs, sOut, kerType, var.hyp, ls_s.hyp, ls_f.hyp, n.starts, n.presample)
      varHyp <- hypers[1]
      lsHyps <- hypers[-1]
    }

    # fill funGpKern slots specific to the functional-input case
    kern@s_lsHyps <- lsHyps[1:ds]
    kern@f_lsHyps <- lsHyps[-c(1:ds)]

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_SF(sMs, fMs, sOut, varHyp, lsHyps[1:ds], lsHyps[(ds+1):(ds+df)], kerType)

    # create objects funGpProj and fill with info specific to the hybrid-input case
    f_proj <- new("funGpProj")
    f_proj@pdims <- f_pdims
    f_proj@family <- f_family
    f_proj@basis <- f_basis
    f_proj@coefs <- f_coefs

    # fill funGp slots specific to the hybrid-input case
    model@ds <- ds
    model@df <- df
    model@f_dims <- f_dims
    model@sIn <- sIn
    model@fIn <- fIn
    model@type = "hybrid"
    model@f_proj <- f_proj

  } else if(!is.null(fIn)) { # functional-input case ***************************************
    # extract information from user inputs specific to the functional-input case
    df <- length(fIn)
    f_dims <- sapply(fIn, ncol)

    # perform projection of functional inputs
    # the projection is such that F = X * B' + e, with
    #
    # n: input points
    # k: original dimension
    # p: projection dimension
    # F: original inputs ............. matrix of dimension nxk
    # B: basis functions ............. matrix of dimension pxk (one basis per column)
    # X: projection coefficients ..... matrix of dimension nxp
    bcj <- dimReduction(fIn, df, f_pdims, f_family)
    f_basis <- bcj$basis
    f_coefs <- bcj$coefs
    f_J <- bcj$J

    # compute functional distance matrices
    fMs <- setFunDistance(f_coefs, f_coefs, f_J)

    # optimize hyperparameters if some is required
    if (all(!is.null(var.hyp), !is.null(ls_f.hyp))) {
      varHyp <- var.hyp
      lsHyps <- ls_f.hyp
    } else {
      hypers <- setHypers_F(fMs, sOut, kerType, var.hyp, ls_f.hyp, n.starts, n.presample)
      varHyp <- hypers[1]
      lsHyps <- hypers[-1]
    }

    # fill funGpKern slots specific to the functional-input case
    kern@f_lsHyps <- lsHyps

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_F(fMs, sOut, varHyp, lsHyps, kerType)

    # create objects funGpProj and fill with info specific to the functional-input case
    f_proj <- new("funGpProj")
    f_proj@pdims <- f_pdims
    f_proj@family <- f_family
    f_proj@basis <- f_basis
    f_proj@coefs <- f_coefs

    # fill funGp slots specific to the functional-input case
    model@ds <- 0
    model@df <- df
    model@f_dims <- f_dims
    model@fIn <- fIn
    model@type = "functional"
    model@f_proj <- f_proj

  } else if(!is.null(sIn)) { # scalar-input case *******************************************
    # extract information from user inputs specific to the scalar-input case
    sIn <- as.matrix(sIn)
    ds <- ncol(sIn)

    # compute scalar distance matrices
    sMs <- setScalDistance(sIn, sIn)

    # optimize hyperparameters if some is required
    if (all(!is.null(var.hyp), !is.null(ls_s.hyp))) {
      varHyp <- var.hyp
      lsHyps <- ls_s.hyp
    } else {
      hypers <- setHypers_S(sIn, sMs, sOut, kerType, var.hyp, ls_s.hyp, n.starts, n.presample)
      varHyp <- hypers[1]
      lsHyps <- hypers[-1]
    }

    # fill funGpKern slots specific to the scalar-input case
    kern@s_lsHyps <- lsHyps

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_S(sMs, sOut, varHyp, lsHyps, kerType)

    # fill funGp slots specific to the scalar-input case
    model@ds <- ds
    model@df <- 0
    model@sIn <- sIn
    model@type = "scalar"

  } else { # error: no inputs were provided
    stop("The user must provide either a scalar-input matrix, a functional-input list or both of them. None has been detected.")
  }

  # fill general funGpKern slots
  kern@kerType <- kerType
  kern@f_disType <- f_disType
  kern@varHyp <- varHyp

  # fill general funGpModel slots
  model@call <- match.call()
  model@sOut <- sOut
  model@n.tot <- n.tr
  model@n.tr <- n.tr
  model@kern <- kern

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
  # 9.    - fpDims ......... array (df) .......... projection dimension of each functional input
  # 10.   - basis .......... list (df) ........... each element (fDims_i x fpDims_i) contains the basis
  #                                                functions used for the projection of one fun. input
  # 11.   - coefs .......... list (df) ........... each element (n x fpDims_i) contains the coefficients
  #                                                used for the projection of one fun. input
  # 12. * kern ............. kernel .............. structures related to the kernel
  # 13.   - kerType ........ char ................ kernel type from {"gauss", "matern5_2", "matern3_2"}
  # 14.   - disType ........ char ................ distance type from {"scalar", "functional"}
  # 15.   - varHyp ......... scalar .............. estimated variance parameter
  # 16.   - lsHyps ......... array (ds + df) ..... estimated length-scale parameters
  # 17. * preMats .......... list (2) ............ KttInv and KttInv.sOut matrices for prediction
  # =====================================================================================================
  return(model)
}
# -------------------------------------------------------------------------------------------------------------------------------------


# Method to print a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name show
#' @description This is my description
#' @rdname show-methods
#' @importFrom methods show
#' @importFrom knitr kable
#' @param object An object to show.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
if(!isGeneric("show")) {setGeneric(name = "show", def = function(object) standardGeneric("show"))}

#' @title Fill!!!!!!!!!!!
#' @name show
#' @rdname show-methods
#' @aliases show,funGp-method
setMethod("show", "funGp", function(object) show.funGp(model = object))

show.funGp <- function(model) {
  mainTxt <- "Gaussian Process Model"
  cat(paste("\n", mainTxt, paste(rep("_", 36), collapse = ""), sep = ""))

  cat(paste("\n\n* Scalar inputs: ", model@ds, "\n", sep = ""))
  cat(paste("* Functional inputs: ", model@df, "", sep = ""))
  if (model@df > 0) {
    # browser()
    np <- min(model@df, 8)
    G <- cbind(paste("F", 1:np, sep = ""), model@f_dims, model@f_proj@pdims, model@f_proj@family, model@kern@f_disType)
    colnames(G) <- c("Input", "Orig. dim", "Proj. dim", "Family", "Distance")
    if (np < model@df) {
      G <- rbind(G, rep("...", 5))
    }
    print(kable(G, align = 'c', row.names = F))
  }

  cat(paste("\n* Total data points: ", model@n.tot, "\n", sep = ""))
  cat(paste("* Trained with: ", model@n.tr, "\n\n", sep = ""))

  cat(paste("* Kernel type: ", model@kern@kerType, "\n", sep = ""))
  cat("* Hyperparameters:\n")
  cat(paste("  -> variance: ", format(model@kern@varHyp, digits = 3, nsmall = 4), "\n", sep = ""))
  cat("  -> length-scale:\n")
  if (model@ds > 0) {
    for (i in 1:model@ds) {
      cat(paste("\t ls(X", i, "): ", format(model@kern@s_lsHyps[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
  }
  if (model@df > 0) {
    for (i in 1:model@df) {
      cat(paste("\t ls(F", i, "): ", format(model@kern@f_lsHyps[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
  }
  cat(paste(rep("_", 58), collapse = ""))
}
# ----------------------------------------------------------------------------------------------------------

# Method to make predictions with a funGp model
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
            predict.funGp(model = object, sIn.pr = sIn.pr, fIn.pr = fIn.pr, detail = detail)
          })

predict.funGp <- function(model, sIn.pr, fIn.pr, detail = "light") {
  # check validity of user inputs
  checkVal_pred_and_sim(as.list(environment()))

  if (all(model@ds > 0, model@df > 0)) { # Hybrid-input case *******************************************
    print("I'm hybrid!")

    # set required data format
    sIn.pr <- as.matrix(sIn.pr)

    # project functional inputs
    f_basis <- model@f_proj@basis
    f_coefs.pr <- mapply(function(B, f) t(solve(crossprod(B), tcrossprod(t(B),f))), f_basis, fIn.pr, SIMPLIFY = F)
    f_J <- lapply(f_basis, crossprod)

    # compute scalar distance matrices
    sMs.tp <- setScalDistance(model@sIn, sIn.pr)
    sMs.pp <- setScalDistance(sIn.pr, sIn.pr)

    # compute functional distance matrices
    fMs.tp <- setFunDistance(model@f_proj@coefs, f_coefs.pr, f_J)
    fMs.pp <- setFunDistance(f_coefs.pr, f_coefs.pr, f_J)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_SF(sMs.tp, sMs.pp, fMs.tp, fMs.pp,
                          model@kern@varHyp, model@kern@s_lsHyps, model@kern@f_lsHyps,
                          model@kern@kerType, model@preMats$L, model@preMats$LInvY, detail)

  } else if (model@df > 0) { # functional-input case *******************************************
    print("I'm functional!")

    # project functional inputs
    f_basis <- model@f_proj@basis
    f_coefs.pr <- mapply(function(B, f) t(solve(crossprod(B), tcrossprod(t(B),f))), f_basis, fIn.pr, SIMPLIFY = F)
    f_J <- lapply(f_basis, crossprod)

    # compute functional distance matrices
    fMs.tp <- setFunDistance(model@f_proj@coefs, f_coefs.pr, f_J)
    fMs.pp <- setFunDistance(f_coefs.pr, f_coefs.pr, f_J)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_F(fMs.tp, fMs.pp, model@kern@varHyp, model@kern@f_lsHyps, model@kern@kerType,
                         model@preMats$L, model@preMats$LInvY, detail)

  } else { # scalar-input case *******************************************
    print("I'm scalar!")

    # set required data format
    sIn.pr <- as.matrix(sIn.pr)

    # compute scalar distance matrices
    sMs.tp <- setScalDistance(model@sIn, sIn.pr)
    sMs.pp <- setScalDistance(sIn.pr, sIn.pr)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_S(sMs.tp, sMs.pp, model@kern@varHyp, model@kern@s_lsHyps, model@kern@kerType,
                         model@preMats$L, model@preMats$LInvY, detail)
  }

  # compute confidence intervals
  preds$lower95 <- preds$mean - qnorm(0.975) * preds$sd
  preds$upper95 <- preds$mean + qnorm(0.975) * preds$sd

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
  return(preds)
}
# ----------------------------------------------------------------------------------------------------------


# Method to make simulations from a funGp model
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
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), matrix(runif(n.tr*22), ncol = 22))
#'
#' # generating output data for training
#' sOut <- as.matrix(sapply(t(1:n.tr), function(i){
#'   x1 <- sIn[i,1]
#'   x2 <- sIn[i,2]
#'   f1 <- fIn[[1]][i,]
#'   f2 <- fIn[[2]][i,]
#'   t1 <- seq(0,1,length = length(f1))
#'   t2 <- seq(0,1,length = length(f2))
#'   as.numeric(x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2))
#' }))
#'
#' # creating a funGp model
#' m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # generating input data for simulation
#' set.seed(100)
#' n.sm <- 100
#' sIn.sm <- expand.grid(x1 = seq(0,1,length = sqrt(n.sm)), x2 = seq(0,1,length = sqrt(n.sm)))
#' fIn.sm <- list(f1 = matrix(runif(n.sm*10), ncol = 10), matrix(runif(n.sm*22), ncol = 22))
#'
#' # making light simulations
#' m1.sims_l <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm)
#'
#' # plotting light simulations
#' plotSims(m1, m1.sims_l)
#'
#' # making full simulations
#' m1.sims_f <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm, detail = "full")
#'
#' # plotting full simulations in full mode
#' plotSims(m1, m1.sims_f)
#'
#' # plotting full simulations in light mode
#' plotSims(m1, m1.sims_f, detail = "light")
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
            simulate.funGp(model = object, nsim = nsim, seed = seed, sIn.sm = sIn.sm, fIn.sm = fIn.sm,
                           nug.sim = nug.sim, detail = detail)
          })

simulate.funGp <- function(model, nsim, seed, sIn.sm, fIn.sm, nug.sim, detail) {
  # check validity of user inputs
  checkVal_pred_and_sim(as.list(environment()))

  # check which type of model it is
  if (all(model@ds > 0, model@df > 0)) { # Hybrid-input case *******************************************
    print("I'm hybrid!")

    # set required data format
    sIn.sm <- as.matrix(sIn.sm)

    # project functional inputs
    f_basis <- model@f_proj@basis
    f_coefs.sm <- mapply(function(B, f) t(solve(crossprod(B), tcrossprod(t(B),f))), f_basis, fIn.sm, SIMPLIFY = F)
    f_J <- lapply(f_basis, crossprod)

    # compute scalar distance matrices
    sMs.ts <- setScalDistance(model@sIn, sIn.sm)
    sMs.ss <- setScalDistance(sIn.sm, sIn.sm)

    # compute functional distance matrices
    fMs.ts <- setFunDistance(model@f_proj@coefs, f_coefs.sm, f_J)
    fMs.ss <- setFunDistance(f_coefs.sm, f_coefs.sm, f_J)

    # make simulations based on the Gaussian Conditioning Theorem
    sims <- makeSims_SF(sMs.ts, sMs.ss, fMs.ts, fMs.ss,
                        model@kern@varHyp, model@kern@s_lsHyps, model@kern@f_lsHyps,
                        model@kern@kerType, model@preMats$L, model@preMats$LInvY, nsim, nug.sim, detail)

  } else if (model@df > 0) { # functional-input case *******************************************
    print("I'm functional!")

    # project functional inputs
    f_basis <- model@f_proj@basis
    f_coefs.sm <- mapply(function(B, f) t(solve(crossprod(B), tcrossprod(t(B),f))), f_basis, fIn.sm, SIMPLIFY = F)
    f_J <- lapply(f_basis, crossprod)

    # compute functional distance matrices
    fMs.ts <- setFunDistance(model@f_proj@coefs, f_coefs.sm, f_J)
    fMs.ss <- setFunDistance(f_coefs.sm, f_coefs.sm, f_J)

    # make simulations based on the Gaussian Conditioning Theorem
    sims <- makeSims_F(fMs.ts, fMs.ss, model@kern@varHyp, model@kern@f_lsHyps, model@kern@kerType,
                       model@preMats$L, model@preMats$LInvY, nsim, nug.sim, detail)

  } else { # scalar-input case *******************************************
    print("I'm scalar!")

    # set required data format
    sIn.sm <- as.matrix(sIn.sm)

    # compute scalar distance matrices
    sMs.ts <- setScalDistance(model@sIn, sIn.sm)
    sMs.ss <- setScalDistance(sIn.sm, sIn.sm)

    # make simulations based on the Gaussian Conditioning Theorem
    sims <- makeSims_S(sMs.ts, sMs.ss, model@kern@varHyp, model@kern@s_lsHyps, model@kern@kerType,
                       model@preMats$L, model@preMats$LInvY, nsim, nug.sim, detail)
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


# Method to update a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name update
#' @rdname update-methods
#' @importFrom stats update
#' @importFrom utils tail
#' @param object An object to update.
#' @param sIn.nw Fill!!!
#' @param fIn.nw Fill!!!
#' @param sOut.nw Fill!!!
#' @param sIn.sb Fill!!!
#' @param fIn.sb Fill!!!
#' @param sOut.sb Fill!!!
#' @param ind.sb Fill!!!
#' @param ind.dl Fill!!!
#' @param var.sb Fill!!!
#' @param ls_s.sb Fill!!!
#' @param ls_f.sb Fill!!!
#' @param var.re Fill!!!
#' @param ls_s.re Fill!!!
#' @param ls_f.re Fill!!!
#' @param ... Further arguments for methods.
#'
# @examples
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export update
setGeneric(name = "update", def = function(object, ...) standardGeneric("update"))

#' @title Updating method for the funGp Class
#' @name update
#' @rdname update-methods
#' @aliases update,funGp-method
setMethod("update", "funGp",
          function(object, sIn.nw = NULL, fIn.nw = NULL, sOut.nw = NULL,
                   sIn.sb = NULL, fIn.sb = NULL, sOut.sb = NULL, ind.sb = NULL,
                   ind.dl = NULL, var.sb = NULL, ls_s.sb = NULL, ls_f.sb = NULL,
                   var.re = F, ls_s.re = F, ls_f.re = F, ...) {
            update.funGp(model = object, sIn.nw = sIn.nw, fIn.nw = fIn.nw, sOut.nw = sOut.nw,
                         sIn.sb = sIn.sb, fIn.sb = fIn.sb, sOut.sb = sOut.sb, ind.sb = ind.sb,
                         ind.dl = ind.dl,
                         var.sb = var.sb, ls_s.sb = ls_s.sb, ls_f.sb = ls_f.sb,
                         var.re = var.re, ls_s.re = ls_s.re, ls_f.re = ls_f.re)
            })

update.funGp <- function(model, sIn.nw, fIn.nw, sOut.nw, sIn.sb, fIn.sb, sOut.sb, ind.sb, ind.dl,
                         var.sb, ls_s.sb, ls_f.sb, var.re, ls_s.re, ls_f.re) {
  # check what does the user want to do
  delInOut <- !is.null(ind.dl)
  subHypers <- any(!is.null(var.sb), !is.null(ls_s.sb), !is.null(ls_f.sb))
  reeHypers <- any(isTRUE(var.re), isTRUE(ls_s.re), isTRUE(ls_f.re))
  if (model@type == "hybrid") {
    subInOut <- any(!is.null(sIn.sb), !is.null(fIn.sb), !is.null(sOut.sb))
    newInOut <- any(!is.null(sIn.nw), !is.null(fIn.nw), !is.null(sOut.nw))
  } else if (model@type == "functional") {
    subInOut <- any(!is.null(fIn.sb), !is.null(sOut.sb))
    newInOut <- any(!is.null(fIn.nw), !is.null(sOut.nw))
  } else if (model@type == "scalar") {
    subInOut <- any(!is.null(sIn.sb), !is.null(sOut.sb))
    newInOut <- any(!is.null(sIn.nw), !is.null(sOut.nw))
  }

  # task names
  # (1) data deletion, (2) data substitution, (3) data addition,
  # (4) var substitution, (5) ls_s substitution, (6) ls_f substitution,
  # (7) var re-estimation, (8) ls_s re-estimation, (9) ls_f re-estimation
  tasknames <- c("data deletion", "data substitution", "data addition",
                 "var substitution", "scalar length-scale substitution", "functional length-scale substitution",
                 "var re-estimation", "scalar length-scale re-estimation", "functional length-scale re-estimation")

  # identify and drop conflicting tasks
  # ----------------------------------------------------
  dptasks <- c()
  if (all(delInOut, subInOut)) { # were deletion and substitution of data both requested?
    dptasks <- c(dptasks, 1, 2)

    if (isTRUE(var.re)) { # was re-estimation of var also requested?
      dptasks <- c(dptasks, 7)
    }
    if (isTRUE(ls_s.re)) { # was re-estimation of ls_s also requested?
      dptasks <- c(dptasks, 8)
    }
    if (isTRUE(ls_f.re)) { # was re-estimation of ls_f also requested?
      dptasks <- c(dptasks, 9)
    }

  }

  if (all(!is.null(var.sb), isTRUE(var.re))) { # were substitution and re-estimation of var both requested?
    dptasks <- c(dptasks, 4, 7)
    var.sb <- NULL
    var.re <- F
  }
  if (all(!is.null(ls_s.sb), isTRUE(ls_s.re))) { # were substitution and re-estimation of ls_s both requested?
    dptasks <- c(dptasks, 5, 8)
    ls_s.sb <- NULL
    ls_s.re <- F
  }
  if (all(!is.null(ls_f.sb), isTRUE(ls_f.re))) { # were substitution and re-estimation of ls_f both requested?
    dptasks <- c(dptasks, 6, 9)
    ls_f.sb <- NULL
    ls_f.re <- F
  }

  if (model@type == "functional") {
    if (isTRUE(ls_s.re)) { # was re-estimation of scalar length-scale coefs requested for a functional model?
      dptasks <- c(dptasks, 8)
    }
  }
  if (model@type == "scalar") {
    if (isTRUE(ls_f.re)) { # was re-estimation of functional length-scale coefs requested for a scalar model?
      dptasks <- c(dptasks, 9)
    }
  }

  # remove duplicates from dropped vector
  dptasks <- unique(dptasks)
  # ----------------------------------------------------

  # perform not dropped tasks
  # ----------------------------------------------------
  modelup <- model
  cptasks <- c()
  if (delInOut & !(1 %in% dptasks)) {
    modelup <- upd_del(model = modelup, ind.dl = ind.dl, remake = all(!newInOut, !subHypers, remake = !reeHypers))
    modelup@call <- model@call
    modelup@n.tr <- model@n.tr
    cptasks <- c(cptasks, 1)
  }
  if (subInOut & !(2 %in% dptasks)) {
    modelup <- upd_subData(model = modelup, sIn.sb = sIn.sb, fIn.sb = fIn.sb,
                           sOut.sb = tryCatch(as.matrix(sOut.sb), error = function(e) sOut.sb), ind.sb = ind.sb,
                           remake = all(!newInOut, !subHypers, !reeHypers))
    modelup@call <- model@call
    cptasks <- c(cptasks, 2)
  }
  if (newInOut) {
    modelup <- upd_add(model = modelup, sIn.nw = sIn.nw, fIn.nw = fIn.nw, sOut.nw = as.matrix(sOut.nw),
                       remake = all(!subHypers, !reeHypers))
    modelup@call <- model@call
    modelup@n.tr <- model@n.tr
    cptasks <- c(cptasks, 3)
  }
  if (subHypers & any(!(c(4,5,6) %in% dptasks))) {
    modelup <- upd_subHypers(model = modelup, var.sb = var.sb, ls_s.sb = ls_s.sb, ls_f.sb = ls_f.sb)
    modelup@call <- model@call
    modelup@n.tr <- model@n.tr
    if (!is.null(var.sb) & !(4 %in% dptasks)) cptasks <- c(cptasks, 4)
    if (!is.null(ls_s.sb) & !(5 %in% dptasks)) cptasks <- c(cptasks, 5)
    if (!is.null(ls_f.sb) & !(6 %in% dptasks)) cptasks <- c(cptasks, 6)
  }
  if (reeHypers & any(!(c(7,8,9) %in% dptasks))) {
    modelup <- upd_reeHypers(model = modelup, var.re = var.re, ls_s.re = ls_s.re, ls_f.re = ls_f.re)
    modelup@call <- model@call
    if (isTRUE(var.re) & !(7 %in% dptasks)) cptasks <- c(cptasks, 7)
    if (isTRUE(ls_s.re) & !(8 %in% dptasks)) cptasks <- c(cptasks, 8)
    if (isTRUE(ls_f.re) & !(9 %in% dptasks)) cptasks <- c(cptasks, 9)
  }
  # ----------------------------------------------------

  # print update summary
  # ----------------------------------------------------
  if (length(cptasks) > 0) { # list of complete tasks if there is any
    cat("--------------\n")
    cat("Update summary\n")
    cat("--------------\n\n")

    cat("* Complete tasks:\n")
    ct <- tasknames[cptasks]
    for (t in ct) {
      cat(paste("  - ", t, "\n", sep = ""))
    }
  }

  if (length(dptasks) > 0) { # list of dropped tasks if there is any
    if (length(cptasks) == 0) {
      cat("--------------\n")
      cat("Update summary\n")
      cat("--------------\n")
    }

    cat("\n* Dropped tasks:\n")
    dt <- tasknames[dptasks]
    for (t in dt) {
      cat(paste("  - ", t, "\n", sep = ""))
    }
    cat("\n* Recall that:\n")
    cat(" - Data points deletion and substitution are not compatible tasks\n")
    cat(" - Hyperparameters substitution and re-estimation are not compatible tasks\n")
    cat(" - Hyperparameters re-estimation is automatically dropped when data deletion and substitution are both requested\n")
    cat(" - Scalar length-scale coeficients re-estimation is automatically dropped when the model has only functional inputs\n")
    cat(" - Functional length-scale coeficients re-estimation is automatically dropped when the model has only scalar inputs\n")
    cat(" -> Please check ?funGp::update for more details\n")
  }
  # ----------------------------------------------------

  return(modelup)
}



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
#' @param ... fill
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod plotLOO
if(!isGeneric("plotLOO")) {setGeneric("plotLOO", function(object, ...) standardGeneric("plotLOO"))}

#' @title Prediction Method for the apk Class
#' @name plotLOO
#' @rdname plotLOO-methods
#' @aliases plotLOO,funGp-method
#' @param xlim something
#' @param ylim something
setMethod("plotLOO", "funGp", function(object, xlim = NULL, ylim = NULL) {
  plotLOO.funGp(model = object, xlim = xlim, ylim = ylim)
  })

plotLOO.funGp <- function(model, xlim, ylim) {
  y_obs <- model@sOut
  R <- tcrossprod(model@preMats$L)/model@kern@varHyp
  Rinv <- solve(R)
  y_pre <- y_obs - diag(Rinv)^(-1) * Rinv %*% y_obs
  xl <- yl <- yr <- range(c(y_obs, y_pre))
  if (!is.null(xlim)) xl <- xlim
  if (!is.null(ylim)) yl <- ylim
  plot(y_obs, y_pre, xlim = xl, ylim = yl, pch = 21, col = "red", bg = "red",
       main = "Model diagnostic by leave-one-out cross-valitation", xlab = "Observed", ylab = "Predicted")
  lines(yr, yr, col = "blue")
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
#' @param xlim_c also
#' @param ylim_c also
#' @param justCal also
#' @param justLin also
setMethod("plotPreds", "funGp",
          function(object, preds, sOut.pr = NULL, xlim_c = NULL, ylim_c = NULL, justCal = F, justLin = F, ...) {
            plotPreds.funGp(preds = preds, sOut.pr = sOut.pr, xlim_c = xlim_c, ylim_c = ylim_c,
                            justCal = justCal, justLin = justLin)
          })

plotPreds.funGp <- function(preds, sOut.pr, xlim_c, ylim_c, justCal, justLin) { # add the usage of user specified limits!!!!!!!!!!!!
  if (all(!is.null(sOut.pr), !justCal, !justLin)) {
    layout(matrix(2:1, nrow = 2))
    par(mar = c(4.1, 4.1, 2.5, 2.1))
  }

  if (!justCal) {
    # sorted mean and 95% limits and true curve
    y <- sort(preds$mean)
    n.pr <- length(y)
    ll <- (preds$lower95)[order(preds$mean)]
    ul <- (preds$upper95)[order(preds$mean)]

    plot(1, type = "n", xlim = c(1, n.pr), ylim = range(y), main = "Sorted predictions", xlab = "Index", ylab = "Predicted")
    x <- 1:n.pr
    polygon(c(x, rev(x)), c(ul, rev(ll)), col = "grey85", border = NA)
    lines(y, col = "red")
    lines(ll, col = "blue")
    lines(ul, col = "blue")
  }

  if (!is.null(sOut.pr)) {
    if (!justCal) {
    # complement for sorted output plot
    lines(sOut.pr[order(preds$mean)], col = "black")
    legend("topleft", legend = c("True", "Pred. mean", "95% CIs"), col = c("black", "red", "blue"), lty = 1, cex = 0.8)
    }

    if (!justLin) {
      # calibration plot
      y_obs <- sOut.pr
      y_pre <- preds$mean
      xl <- yl <- yr <- range(c(y_obs, y_pre))
      if (!is.null(xlim_c)) xl <- xlim_c
      if (!is.null(ylim_c)) yl <- ylim_c
      plot(y_obs, y_pre, xlim = xl, ylim = yl, pch = 21, col = "red", bg = "red",
           main = "Model predictions at new input points", xlab = "Observed", ylab = "Predicted")
      lines(y_obs, y_obs, col = "blue")
    }
  } else {
    # complement for sorted output plot
    lines(sOut.pr[order(preds$mean)], col = "black")
    legend("topleft", legend = c("Pred. mean", "95% CIs"), col = c("black", "red", "blue"), lty = 1, cex = 0.8)
  }

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
    matplot(t(sims), type = "l", lty = 1, col = "palegreen4",
            main = "Simulations from a funGp model", xlab = "Sim. index", ylab = "Output")
  } else {
    if (detail == "light") {
      matplot(t(sims$obs), type = "l", lty = 1, col = "palegreen4",
              main = "Simulations from a funGp model", xlab = "Sim. index", ylab = "Output")
    } else {
      matplot(t(sims$obs), type = "l", lty = 1, col = "grey",
              main = "Simulations from a funGp model", xlab = "Sim. index", ylab = "Output")
      lines(sims$mean, col = "red")
      lines(sims$lower95, col = "blue")
      lines(sims$upper95, col = "blue")
      legend("topleft", legend = c("Sims", "Mean", "95% CIs"), col = c("grey50", "red", "blue"), lty = 1, cex = 0.8)
    }
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
setMethod("getCoef", "funGp", function(object) getCoef.funGp(model = object))

getCoef.funGp <- function(model) {
  coefs <- model@kern@varHyp
  names_ls_s <- c()
  if (model@ds > 0) {
    names_ls_s <- paste("ls(X", 1:model@ds, ")", sep = "")
    coefs <- c(coefs, model@kern@s_lsHyps)
  }
  names_ls_f <- c()
  if (model@df > 0) {
    names_ls_f <- paste("ls(F", 1:model@df, ")", sep = "")
    coefs <- c(coefs, model@kern@f_lsHyps)
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
setMethod("getBasis", "funGp", function(object) getBasis.funGp(model = object))

getBasis.funGp <- function(model) {
  if (model@df > 0) {
    return(model@proj@basis)
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
setMethod("getRedfIn", "funGp", function(object) getRedfIn.funGp(model = object))

getRedfIn.funGp <- function(model) {
  if (model@df > 0) {
    rfIn <- model@proj@coefs
    names(rfIn) <- paste("alpha(F", 1:model@df, ")", sep = "")
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
setMethod("getProjfIn", "funGp", function(object) getProjfIn.funGp(model = object))

getProjfIn.funGp <- function(model) {
  if (model@df > 0) {
    fpIn <- mapply(function(m1, m2) m1 %*% t(m2), m1 = model@proj@coefs, m2 = model@proj@basis)
    names(fpIn) <- paste("Fp", 1:model@df, sep = "")
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
setMethod("getProjgram", "funGp", function(object) getProjgram.funGp(model = object))

getProjgram.funGp <- function(model) {
  if (model@df > 0) {
    gram <- list()
    for (i in 1:model@df) {
      gram[[i]] <- t(model@proj@basis[[i]]) %*% model@proj@basis[[i]]
    }
    names(gram) <- paste("G(F", 1:model@df, ")", sep = "")
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
setMethod("getTrainCov", "funGp", function(object) getTrainCov.funGp(model = object))

getTrainCov.funGp <- function(model) {
  return(tcrossprod(model@preMats$L))
}
# ----------------------------------------------------------------------------------------------------------
