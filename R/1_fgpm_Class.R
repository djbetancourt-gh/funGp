# ==========================================================================================================
# S4 class for funGp Gaussian process models
# ==========================================================================================================
#' @title S4 class for funGp Gaussian process models
#' @description This is the formal representation of Gaussian process models within the
#'   \link[=funGp-package]{funGp package}. Gaussian process models are useful statistical tools in the
#'   modeling of complex input-output relationships.
#' \itemize{
#'  \item{\strong{Main methods}}{\cr
#'        \link[funGp]{fgpm}: creation of funGp regression models \cr
#'        \link[funGp]{predict,fgpm-method}: output estimation at new input points based on a \code{fgpm} model \cr
#'        \link[funGp]{simulate,fgpm-method}: random sampling from a \code{fgpm} model \cr
#'        \link[funGp]{update,fgpm-method}: modification of data and hyperparameters of a \code{fgpm} model
#'  }
#'  \item{\strong{Plotters}}{\cr
#'         \link[funGp]{plot,fgpm-method}: validation plot for a \code{fgpm} model \cr
#'        \link[funGp]{plot.predict.fgpm}: plot of predictions based on a \code{fgpm} model \cr
#'        \link[funGp]{plot.simulate.fgpm}: plot of simulations based on a \code{fgpm} model
#'  }
#' }
#'
#' @slot howCalled Object of class \code{"\linkS4class{modelCall}"}. User call reminder.
#' @slot type Object of class \code{"character"}. Type of model based on type of inputs. To be set from
#'   {"scalar", "functional", "hybrid"}.
#' @slot ds Object of class \code{"numeric"}. Number of scalar inputs.
#' @slot df Object of class \code{"numeric"}. Number of functional inputs.
#' @slot f_dims Object of class \code{"numeric"}. An array with the original dimension of each functional
#'   input.
#' @slot sIn Object of class \code{"matrix"}. The scalar input points. Variables are arranged by columns and
#'   coordinates by rows.
#' @slot fIn Object of class \code{"list"}. The functional input points. Each element of the list contains
#'   a functional input in the form of a matrix. In each matrix, curves representing functional coordinates
#'   are arranged by rows.
#' @slot sOut Object of class \code{"matrix"}. The scalar output values at the coordinates specified by sIn
#'   and/or fIn.
#' @slot n.tot Object of class \code{"integer"}. Number of observed points used to compute the training-training
#'   and training-prediction covariance matrices.
#' @slot n.tr Object of class \code{"integer"}. Among all the points loaded in the model, the amount used for
#'   training.
#' @slot f_proj Object of class \code{"fgpProj"}. Data structures related to the projection of functional
#'   inputs. Check \linkS4class{fgpProj} for more details.
#' @slot kern Object of class \code{"fgpKern"}. Data structures related to the kernel of the Gaussian process
#'   model. Check \linkS4class{fgpKern} for more details.
#' @slot nugget Object of class \code{"numeric"}. Variance parameter standing for the homogeneous nugget effect.
#' @slot preMats Object of class \code{"list"}. L and LInvY matrices pre-computed for prediction. L is a lower
#'   diagonal matrix such that \eqn{L'L} equals the training auto-covariance matrix \eqn{K.tt}. On the other
#'   hand, \eqn{LInvY = L^(-1) * sOut}.
#' @slot convergence Object of class \code{"numeric"}. Integer code either confirming convergence or indicating
#'  an error. Check the convergence component of the Value returned by \code{\link[stats]{optim}}.
#' @slot NegLogLike Object of class \code{"numeric"}. Negated log-likelihood obained by \code{\link[stats]{optim}}
#'  during hyperparameter optimization.
#'
#' @section Useful material:
#' \itemize{
#'  \item{\strong{Manual}}{
#'  \href{https://hal.archives-ouvertes.fr/hal-02536624}{
#'  Gaussian Process Regression for Scalar and Functional Inputs with funGp - The in-depth tour
#'  }}
#' }
#'
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
#'
#' @include 2_fgpProj_Class.R
#' @include 2_fgpKern_Class.R
#' @include 8_outilsCode.R
#'
#' @rdname fgpm-class
#' @export
setClass("fgpm",
         slots = c(
           howCalled = "modelCall",    # reminder of function call
           type = "character",         # type of model. To be set from {"scalar", "functional", "hybrid"}.
           ds = "numeric",             # number of scalar inputs
           df = "numeric",             # number of functional inputs
           f_dims = "numeric",         # dimension of each functional input
           sIn = "matrix",             # scalar inputs
           fIn = "list",               # each element (n x fDims_i) contains a functional input
           sOut = "matrix",            # scalar output
           n.tot = "integer",          # number of observed points
           n.tr = "integer",           # number of observed points used for training
           f_proj = "fgpProj",         # structures related to the projection of functional inputs
           kern = "fgpKern",           # structures related to the kernel of the model
           nugget = "numeric",         # variance parameter standing for the homogeneous nugget effect
           preMats = "list",           # pre-computed KttInv and KttInv.sOut matrices
           convergence = "numeric",    # indicator of convergence/error
           NegLogLike = "numeric"      # negated loglikelihood achieved during hypers optim
         ),
         validity = function(object) {TRUE})
# ==========================================================================================================



# ==========================================================================================================
# Fitting of a funGp model
# ==========================================================================================================
#' @title Gaussian process models for scalar and functional inputs
#' @description This function enables fitting of Gaussian process regression models. The inputs can be
#'   either scalar, functional or a combination of both types.
#'
#' @param sIn An optional matrix of scalar input values to train the model. Each column must match an input
#'   variable and each row a training point. Either scalar input coordinates (sIn), functional input
#'   coordinates (fIn), or both must be provided.
#' @param fIn An optional list of functional input values to train the model. Each element of the list must
#'   be a matrix containing the set of curves corresponding to one functional input. Either scalar input
#'   coordinates (sIn), functional input coordinates (fIn), or both must be provided.
#' @param sOut A vector (or 1-column matrix) containing the values of the scalar output at the specified
#'   input points.
#' @param kerType An optional character string specifying the covariance structure to be used. To be chosen
#'   between "gauss", "matern5_2" and "matern3_2". Default is "matern5_2".
#' @param f_disType An optional array of character strings specifying the distance function to be used for
#'   each functional coordinates within the covariance function of the Gaussian process. To be chosen between
#'   "L2_bygroup" and "L2_byindex". The L2_bygroup distance considers each curve as a whole and uses a single
#'   length-scale parameter per functional input variable. The L2_byindex distance uses as many length-scale
#'   parameters per functional input as discretization points it has. For instance an input discretized as
#'   a vector of size 8 will use 8 length-scale parameters when using L2_byindex. If dimension reduction of
#'   a functional input is requested, then L2_byindex uses as many length scale parameters as effective
#'   dimensions used to represent the input. A single character string can also be passed as a general
#'   selection for all the functional inputs of the model. More details in
#'   \href{https://www.sciencedirect.com/science/article/abs/pii/S0951832019301693}{
#'   the reference article}
#'   and
#'   \href{https://hal.archives-ouvertes.fr/hal-02536624}{
#'   the in-depth package manual}. Default is "L2_bygroup".
#' @param f_pdims An optional array with the projection dimension for each functional input. For each input,
#'   the projection dimension should be an integer between 0 and its original dimension, with 0 denoting
#'   no projection. A single character string can also be passed as a general selection for all the functional
#'   inputs of the model. Default is 3.
#' @param f_basType An optional array of character strings specifying the family of basis functions to be used
#'   in the projection of each functional input. To be chosen between "B-splines" and "PCA". A single character
#'   string can also be passed as a general selection for all the functional inputs of the model. This argument
#'   will be ignored for those inputs for which no projection was requested (i.e., for which f_pdims = 0).
#'   Default is "B-splines".
#' @param var.hyp An optional number indicating the value that should be used as the variance parameter of the
#'   model. If not provided, it is estimated through likelihood maximization.
#' @param ls_s.hyp An optional numeric array indicating the values that should be used as length-scale parameters
#'   for the scalar inputs. If provided, the size of the array should match the number of scalar inputs. If not
#'   provided, these parameters are estimated through likelihood maximization.
#' @param ls_f.hyp An optional numeric array indicating the values that should be used as length-scale parameters
#'   for the functional inputs. If provided, the size of the array should match the number of effective dimensions.
#'   Each input using the "L2_bygroup" distance will count 1 effective dimension, and each input using the
#'   "L2_byindex" distance will count as many effective dimensions as specified by the corresponding element of
#'   the f_pdims argument. For instance, two functional inputs of original dimensions 10 and 22, the first one
#'   projected onto a space of dimension 5 with "L2_byindex" distance, and the second one not projected with
#'   "L2_bygroup" distance will make up a total of 6 effective dimensions; five for the first functional input and
#'   one for second one. If this argument is not provided, the functional length-scale parameters are estimated
#'   through likelihood maximization.
#' @param nugget An optional variance value standing for the homogeneous nugget effect. A tiny nugget might help
#'   to overcome numerical problems related to the ill-conditioning of the covariance matrix. Default is 1e-8.
#' @param n.starts An optional integer indicating the number of initial points to use for the optimization of the
#'   hyperparameters. A parallel processing cluster can be exploited in order to speed up the evaluation of
#'   multiple initial points. More details in the description of the argument par.clust below. Default is 1.
#' @param n.presample An optional integer indicating the number of points to be tested in order to select the
#'   n.starts initial points. The n.presample points will be randomly sampled from the hyper-rectangle defined by: \cr \cr
#'   1e-10 \eqn{\le} \code{ls_s.hyp[i]} \eqn{\le} 2*max(\code{sMs[[i]]}), for i in 1 to the number of scalar inputs, \cr
#'   1e-10 \eqn{\le} \code{ls_f.hyp[i]} \eqn{\le} 2*max(\code{fMs[[i]]}), for i in 1 to the number of functional inputs, \cr \cr
#'   with  sMs and fMs the lists of distance matrices for the scalar and functional inputs, respectively. The value of
#'   n.starts will be assigned to n.presample if this last is smaller. Default is 20.
#' @param par.clust An optional parallel processing cluster created with the \code{\link[parallel]{makeCluster}} function
#'   of the \link[=parallel]{parallel package}. If not provided, multistart optimizations are done in sequence.
#' @param trace An optional boolean indicating if control messages from the \link[stats]{optim} function regarding the
#'   optimization of the hyperparameters should be printed to console. Default is TRUE.
#' @param pbars An optional boolean indicating if progress bars should be displayed. Default is TRUE.
#' @param control.optim An optional list to be passed as the \code{control} argument to \code{\link[stats]{optim}}, the function
#'   in charge of the non-linear optimization of the hyperparameters. Default is \code{list(trace = TRUE)}, equivalent to
#'   \code{list(trace = 1)}, which enables the printing of tracing information on the progress of the optimization. Before
#'   interacting with this \code{\link[funGp]{fgpm}} argument, please carefully check the documentation provided in
#'   \code{\link[stats]{optim}} to ensure a coherent behavior and sound results. Note that: (i) at this time, only the
#'   \code{"L-BFGS-B"} method (Byrd et. al., 1995) is enabled in \code{\link[funGp]{fgpm}}; (ii) \code{control.optim$fnscale}
#'   should not be used since our optimization problem is strictly of minimization, not maximization.
#'
#' @return An object of class \linkS4class{fgpm} containing the data structures representing the fitted funGp model.
#'
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
#'
#' @references Betancourt, J., Bachoc, F., Klein, T., Idier, D., Pedreros, R., and Rohmer, J. (2020),
#' "Gaussian process metamodeling of functional-input code for coastal flood hazard assessment".
#' \emph{Reliability Engineering & System Safety}, \strong{198}, 106870.
#' \href{https://www.sciencedirect.com/science/article/abs/pii/S0951832019301693}{[RESS]}
#' \href{https://hal.archives-ouvertes.fr/hal-01998724}{[HAL]}
#'
#' @references Betancourt, J., Bachoc, F., Klein, T., and Gamboa, F. (2020),
#' Technical Report: "Ant Colony Based Model Selection for Functional-Input Gaussian Process Regression. Ref. D3.b (WP3.2)".
#' \emph{RISCOPE project}.
#' \href{https://hal.archives-ouvertes.fr/hal-02532713}{[HAL]}
#'
#' @references Betancourt, J., Bachoc, F., and Klein, T. (2020),
#' R Package Manual: "Gaussian Process Regression for Scalar and Functional Inputs with funGp - The in-depth tour".
#' \emph{RISCOPE project}.
#' \href{https://hal.archives-ouvertes.fr/hal-02536624}{[HAL]}
#'
#' @seealso \strong{*} \link[funGp]{plot,fgpm-method}: validation plot for a \code{fgpm} model;
#' @seealso \strong{*} \link[funGp]{predict,fgpm-method} for predictions based on a \code{fgpm} model;
#' @seealso \strong{*} \link[funGp]{simulate,fgpm-method} for simulations based on a \code{fgpm} model;
#' @seealso \strong{*} \link[funGp]{update,fgpm-method} for post-creation updates on a \code{fgpm} model;
#' @seealso \strong{*} \link[funGp]{fgpm_factory} for funGp heuristic model selection.
#'
#' @examples
#' # creating funGp model using default fgpm arguments________________________________________
#' # generating input data for training
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#'
#' # generating output data for training
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#'
#' # building a scalar-input funGp model
#' ms <- fgpm(sIn = sIn, sOut = sOut)
#'
#' # building a functional-input funGp model
#' mf <- fgpm(fIn = fIn, sOut = sOut)
#'
#' # building a hybrid-input funGp model
#' msf <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # plotting the three models
#' plot(ms)
#' plot(mf)
#' plot(msf)
#'
#' # printing the three models
#' summary(ms) # equivalent to show(ms)
#' summary(mf) # equivalent to show(mf)
#' summary(msf) # equivalent to show(msf)
#'
#'
#' # recovering useful information from a funGp model_________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # recovering data from model slots
#' m1@f_proj@coefs # list of projection coefficients for the functional inputs
#' m1@f_proj@basis # list of projection basis functions for the functional inputs
#' Map(function(a, b) a %*% t(b), m1@f_proj@coefs, m1@f_proj@basis) # list of projected
#'                                                                  # functional inputs
#' tcrossprod(m1@preMats$L) # training auto-covariance matrix
#'
#'
#' # making predictions based on a funGp model________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # generating input data for prediction
#' n.pr <- 100
#' sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
#'                                 x2 = seq(0,1,length = sqrt(n.pr))))
#' fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), matrix(runif(n.pr*22), ncol = 22))
#'
#' # making predictions
#' m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
#'
#' # plotting predictions
#' plot(m1.preds)
#'
#'
#' # simulating from a funGp model____________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # generating input data for simulation
#' n.sm <- 100
#' sIn.sm <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.sm)),
#'                                 x2 = seq(0,1,length = sqrt(n.sm))))
#' fIn.sm <- list(f1 = matrix(runif(n.sm*10), ncol = 10), matrix(runif(n.sm*22), ncol = 22))
#'
#' # making simulations
#' m1.sims <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm)
#'
#' # plotting simulations
#' plot(m1.sims)
#'
#'
#' # creating funGp model using custom fgpm arguments_________________________________________
#' # generating input and output data
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#'
#' # original dimensions
#' # f1: 10
#' # f2: 22
#'
#' # building a the model with the following structure
#' #    - Kernel: Gaussian
#' #    - f1: L2_byindex distance, no projection -> 10 length-scale parameters
#' #    - f2: L2_bygroup distance, B-spline basis of dimension 5 -> 1 length-scale parameter
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut,
#'            kerType = "gauss", f_disType = c("L2_byindex", "L2_bygroup"),
#'            f_pdims = c(0,5), f_basType = c(NA, "B-splines"))
#'
#' # plotting the model
#' plot(m1)
#'
#' # printing the model
#' m1 # equivalent to show(m1)
#'
#' \dontrun{
#' # multistart and parallelization in fgpm___________________________________________________
#' # generating input and output data
#' set.seed(100)
#' n.tr <- 243
#' sIn <- expand.grid(x1 = seq(0,1,length = n.tr^(1/5)), x2 = seq(0,1,length = n.tr^(1/5)),
#'                    x3 = seq(0,1,length = n.tr^(1/5)), x4 = seq(0,1,length = n.tr^(1/5)),
#'                    x5 = seq(0,1,length = n.tr^(1/5)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB7(sIn, fIn, n.tr)
#'
#' # calling fgpm with multistart in parallel
#' cl <- parallel::makeCluster(2)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut, n.starts = 10, par.clust = cl) # (~14 seconds)
#' parallel::stopCluster(cl)
#'
#' # NOTE: in order to provide progress bars for the monitoring of time consuming processes
#' #       ran in parallel, funGp relies on the doFuture and future packages. Parallel processes
#' #       suddenly interrupted by the user tend to leave corrupt connections. This problem is
#' #       originated outside funGp, which limits our control over it. In the manual
#' #       of funGp, we provide a temporary solution to the issue and we remain attentive in
#' #       case it appears a more elegant way to handle it or a manner to suppress it.
#' #
#' #       funGp manual: https://hal.archives-ouvertes.fr/hal-02536624
#' }
#'
#' @importFrom methods new
#' @importFrom microbenchmark microbenchmark
#' @export
fgpm <- function(sIn = NULL, fIn = NULL, sOut, kerType = "matern5_2",
                 f_disType = "L2_bygroup", f_pdims = 3, f_basType = "B-splines",
                 var.hyp = NULL, ls_s.hyp = NULL, ls_f.hyp = NULL, nugget = 1e-8,
                 n.starts = 1, n.presample = 20, par.clust = NULL, trace = TRUE, pbars = TRUE,
                 control.optim = list(trace = TRUE)) {
  # extend simplified user inputs to full versions
  if (!is.null(sIn)) {
    if (is.numeric(sIn)) sIn <- as.matrix(sIn)
  }
  if (!is.null(fIn)) {
    if (is.matrix(fIn)) fIn <- list(fIn)
    if (length(f_disType) == 1) f_disType <- rep(f_disType, length(fIn))
    if (length(f_pdims) == 1) f_pdims <- rep(f_pdims, length(fIn))
    if (length(f_basType) == 1) f_basType <- rep(f_basType, length(fIn))
    f_basType[which(f_pdims == 0)] <- NA
  }

  # check validity of user inputs
  checkVal_fgpm(as.list(environment()))

  # create objects of class fgpKern and fgpm
  kern <- new("fgpKern")
  model <- new("fgpm")

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
    # n: input points
    # k: original dimension
    # p: projection dimension
    # F: original inputs ............. matrix of dimension nxk
    # B: basis functions ............. matrix of dimension pxk (one basis per column)
    # X: projection coefficients ..... matrix of dimension nxp
    bcj <- dimReduction(fIn, df, f_pdims, f_basType)
    f_basis <- bcj$basis
    f_coefs <- bcj$coefs
    f_J <- bcj$J

    # compute scalar distance matrices
    sMs <- setDistMatrix_S(sIn, sIn)

    # compute functional distance matrices
    fMs <- setDistMatrix_F(f_coefs, f_coefs, f_J, f_disType)
    owners <- paste("F", getOwners(df, f_disType, sapply(f_coefs, ncol)), sep = "")

    # optimize hyperparameters if some is required
    if (all(!is.null(var.hyp), !is.null(ls_s.hyp), !is.null(ls_f.hyp))) {
      varHyp <- var.hyp
      lsHyps <- c(ls_s.hyp, ls_f.hyp)
    } else {
      optResult <- setHypers_SF(sMs, fMs, sOut, kerType, var.hyp, ls_s.hyp, ls_f.hyp,
                                n.starts, n.presample, nugget, par.clust, trace, pbars,
                                control.optim)
      hypers <- optResult$hypers
      varHyp <- hypers[1]
      lsHyps <- hypers[-1]
    }

    # fill fgpKern slots specific to the functional-input case
    kern@s_lsHyps <- lsHyps[1:ds]
    kern@f_lsHyps <- lsHyps[-c(1:ds)]
    kern@f_lsOwners <- owners

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_SF(sMs, fMs, sOut, varHyp, lsHyps[1:ds], lsHyps[-c(1:ds)], kerType, nugget)

    # create objects fgpProj and fill with info specific to the hybrid-input case
    f_proj <- new("fgpProj")
    f_proj@pdims <- f_pdims
    f_proj@basType <- f_basType
    f_proj@basis <- f_basis
    f_proj@coefs <- f_coefs

    # fill fgpm slots specific to the hybrid-input case
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
    bcj <- dimReduction(fIn, df, f_pdims, f_basType)
    f_basis <- bcj$basis
    f_coefs <- bcj$coefs
    f_J <- bcj$J

    # compute functional distance matrices
    fMs <- setDistMatrix_F(f_coefs, f_coefs, f_J, f_disType)
    owners <- paste("F", getOwners(df, f_disType, sapply(f_coefs, ncol)), sep = "")

    # optimize hyperparameters if some is required
    if (all(!is.null(var.hyp), !is.null(ls_f.hyp))) {
      varHyp <- var.hyp
      lsHyps <- ls_f.hyp
    } else {
      optResult <- setHypers_F(fMs, sOut, kerType, var.hyp, ls_f.hyp,
                               n.starts, n.presample, nugget, par.clust, trace, pbars,
                               control.optim)
      hypers <- optResult$hypers
      varHyp <- hypers[1]
      lsHyps <- hypers[-1]
    }

    # fill fgpKern slots specific to the functional-input case
    kern@f_lsHyps <- lsHyps
    kern@f_lsOwners <- owners

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_F(fMs, sOut, varHyp, lsHyps, kerType, nugget)

    # create objects fgpProj and fill with info specific to the functional-input case
    f_proj <- new("fgpProj")
    f_proj@pdims <- f_pdims
    f_proj@basType <- f_basType
    f_proj@basis <- f_basis
    f_proj@coefs <- f_coefs

    # fill fgpm slots specific to the functional-input case
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
    sMs <- setDistMatrix_S(sIn, sIn)

    # optimize hyperparameters if some is required
    if (all(!is.null(var.hyp), !is.null(ls_s.hyp))) {
      varHyp <- var.hyp
      lsHyps <- ls_s.hyp
    } else {
      optResult <- setHypers_S(sIn, sMs, sOut, kerType, var.hyp, ls_s.hyp,
                               n.starts, n.presample, nugget, par.clust, trace, pbars,
                               control.optim)
      hypers <- optResult$hypers
      varHyp <- hypers[1]
      lsHyps <- hypers[-1]
    }

    # fill fgpKern slots specific to the scalar-input case
    kern@s_lsHyps <- lsHyps

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_S(sMs, sOut, varHyp, lsHyps, kerType, nugget)

    # fill fgpm slots specific to the scalar-input case
    model@ds <- ds
    model@df <- 0
    model@sIn <- sIn
    model@type = "scalar"

  } else { # error: no inputs were provided
    stop("The user must provide either a scalar-input matrix, a functional-input list or both of them. None has been detected.")
  }

  # fill general fgpKern slots
  kern@kerType <- kerType
  kern@f_disType <- f_disType
  kern@varHyp <- varHyp

  # fill general fgpm slots
  model@howCalled@string <- gsub("^ *|(?<= ) | *$", "", paste0(deparse(match.call()), collapse = " "), perl = TRUE)
  model@sOut <- sOut
  model@n.tot <- n.tr
  model@n.tr <- n.tr
  model@kern <- kern
  model@nugget <- nugget
  if (exists("optResult")) {
    model@convergence <- optResult$convg
    model@NegLogLike <- optResult$nllik
  } else {
    model@convergence <- as.numeric(NA)
    model@NegLogLike <- as.numeric(NA)
  }

  # ________________________________________________________________________________________________________
  # Attributes checklist
  # ________________________________________________________________________________________________________
  # 0.  * howCalled ........ call ................ functional call
  # 1.  * type ............. char ................ type of inputs from {"scalar", "functional", "hybrid"}
  # 2.  * ds ............... scalar .............. number of scalar inputs
  # 3.  * df ............... scalar .............. number of functional inputs
  # 4.  * fDims ............ array (df) .......... dimension of each functional input
  # 5.  * sIn .............. matrix (n x ds) ..... scalar inputs
  # 6.  * fIn .............. list (df) ........... each element (n x fDims_i) contains a functional input
  # 7.  * sOut ............. matrix (n x 1) ...... scalar output
  # 8.  * n.tot ............ scalar .............. total number of points loaded in the model
  # 9.  * n.tr ............. scalar .............. among all the loaded points, the amount used for training
  # 10. * f_proj ........... proj ................ structures related to the projection of fun. inputs
  # 11.   - pdims .......... array (df) .......... projection dimension of each functional input
  # 12.   - basType ........ array (df) .......... family of basis functions used for each input
  # 13.   - basis .......... list (df) ........... each element (fDims_i x fpDims_i) contains the basis
  #                                                functions used for the projection of one fun. input
  # 14.   - coefs .......... list (df) ........... each element (n x fpDims_i) contains the coefficients
  #                                                used for the projection of one fun. input
  # 15. * kern ............. kernel .............. structures related to the kernel
  # 16.   - kerType ........ char ................ kernel type from {"gauss", "matern5_2", "matern3_2"}
  # 17.   - f_disType ...... char ................ distance type from {"scalar", "functional"}
  # 18.   - varHyp ......... scalar .............. estimated variance parameter
  # 19.   - s_lsHyps ....... array (ds) .......... estimated length-scale parameters for scl. inputs
  # 20.   - f_lsHyps ....... array (ds) .......... estimated length-scale parameters for fun. inputs
  # 21.   - f_lsOwners ..... array (ds) .......... input linked to each fun length-scale coefficient
  # 22. * nugget ........... scalar .............. homogeneous nugget effect
  # 23. * preMats .......... list (2) ............ KttInv and KttInv.sOut matrices for prediction
  # ________________________________________________________________________________________________________
  return(model)
}
# ==========================================================================================================



# ==========================================================================================================
# Printing of a \code{fgpm} model
# ==========================================================================================================
#' @importFrom knitr kable
#' @rdname show-methods
#' @aliases show,fgpm-method
#' @noRd
setMethod("show", "fgpm", function(object) show.fgpm(model = object))

show.fgpm <- function(model) {
  mainTxt <- "Gaussian Process Model"
  if (model@df > 0) {
    cat(paste("\n", mainTxt, paste(rep("_", 36), collapse = ""), sep = ""))
  } else {
    cat(paste("\n", mainTxt, paste(rep("_", 2), collapse = ""), sep = ""))
  }

  cat(paste("\n* Scalar inputs: ", model@ds, "\n", sep = ""))
  cat(paste("* Functional inputs: ", model@df, "", sep = ""))
  if (model@df > 0) {
    np <- min(model@df, 8)
    G <- cbind(paste("F", 1:np, sep = ""), model@f_dims, model@f_proj@pdims, model@f_proj@basType, model@kern@f_disType)
    colnames(G) <- c("Input", "Orig. dim", "Proj. dim", "Basis", "Distance")
    if (np < model@df) {
      G <- rbind(G, rep("...", 5))
    }
    print(kable(G, align = 'c', row.names = FALSE))
  }

  cat(paste("\n* Total data points: ", model@n.tot, "\n", sep = ""))
  cat(paste("* Trained with: ", model@n.tr, "\n", sep = ""))

  cat(paste("* Kernel type: ", model@kern@kerType, "\n", sep = ""))
  cat(paste("* Convergence: ", model@convergence, "\n", sep = ""))
  cat(paste("* NegLogLik: ", format(model@NegLogLike, digits = 3, nsmall = 4), "\n", sep = ""))
  cat("* Hyperparameters:\n")
  cat(paste("  -> variance: ", format(model@kern@varHyp, digits = 3, nsmall = 4), "\n", sep = ""))
  cat("  -> length-scale:\n")
  max.pr <- 8
  if (model@type == "hybrid") {
    # prepare lenght-scale parameters for printing (allows to print a maximum of 8 length-scale parameters)
    a <- paste("X", 1:model@ds, sep = "")
    b <- model@kern@f_lsOwners
    all.owners <- c(a, b)[order(c(seq_along(a)*2 - 1, seq_along(b)*2))]
    a <- model@kern@s_lsHyps
    b <- model@kern@f_lsHyps
    all.ls <- c(a, b)[order(c(seq_along(a)*2 - 1, seq_along(b)*2))]
    top.owners <- all.owners[1:min(max.pr,length(all.owners))]
    top.ls <- all.ls[1:min(max.pr,length(all.owners))]
    ids.s <- grepl("X", top.owners)
    ids.f <- grepl("F", top.owners)
    s.ls <- top.ls[ids.s]
    s.owners <- top.owners[ids.s]
    f.ls <- top.ls[ids.f]
    f.owners <- top.owners[ids.f]

    for (i in 1:length(s.ls)) {
      cat(paste("\t ls(", s.owners[i], "): ", format(s.ls[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
    for (i in 1:length(f.ls)) {
      cat(paste("\t ls(", f.owners[i], "): ", format(f.ls[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
    if (length(all.ls) > max.pr)
      cat("\n Some length-scale parameters were not printed. Consider\n checking 'model@kern@s_lsHyps' and 'model@kern@f_lsHyps'\n")

  } else if (model@type == "functional") {
    ids.f <- 1:min(max.pr,length(model@kern@f_lsHyps))
    f.ls <- model@kern@f_lsHyps[ids.f]
    f.owners <- model@kern@f_lsOwners[ids.f]
    for (i in 1:length(f.ls)) {
      cat(paste("\t ls(", f.owners[i], "): ", format(f.ls[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
    if (length(f.ls) > max.pr)
      cat("\n Some length-scale parameters were not printed. Consider\n checking 'model@kern@f_lsHyps'\n")

  } else {
    ids.s <- 1:min(max.pr,length(model@kern@s_lsHyps))
    s.ls <- model@kern@s_lsHyps[ids.s]
    s.owners <- paste("X", 1:model@ds, sep = "")
    for (i in 1:length(s.ls)) {
      cat(paste("\t ls(", s.owners[i], "): ", format(s.ls[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
    if (length(s.ls) > max.pr)
      cat("\n Some length-scale parameters were not printed. Consider\n checking 'model@kern@s_lsHyps'\n")
  }
  if (model@df > 0) {
    cat(paste(rep("_", 58), collapse = ""), "\n")
  } else {
    cat(paste(rep("_", 24), collapse = ""), "\n")
  }
}
# ==========================================================================================================



# ==========================================================================================================
# Prediction based on a \code{fgpm} model
# ==========================================================================================================
#' @name predict
#' @rdname predict-methods
#' @importFrom stats predict
#' @exportMethod predict
#' @noRd
setGeneric(name = "predict", def = function(object, ...) standardGeneric("predict"))

#' @title Prediction from a \code{fgpm} Gaussian process model
#' @description This method enables prediction based on a \code{fgpm} model, at any given set of
#'   points. Check \code{\link{fgpm}} for information on how to create \code{fgpm} models.
#'
#' @param object An object of class \linkS4class{fgpm} corresponding to the funGp model that should be used
#'   to predict the output.
#' @param ... Not used.
#' @param sIn.pr An optional matrix of scalar input coordinates at which the output values should be
#'   predicted. Each column is interpreted as a scalar input variable and each row as a coordinate.
#'   Either scalar input coordinates (sIn.pr), functional input coordinates (fIn.pr), or both must be provided.
#' @param fIn.pr An optional list of functional input coordinates at which the output values should be
#'   predicted. Each element of the list is interpreted as a functional input variable. Every functional input
#'   variable should be provided as a matrix with one curve per row. Either scalar input coordinates (sIn.pr),
#'   functional input coordinates (fIn.pr), or both must be provided.
#' @param detail An optional character specifying the extent of information that should be delivered
#'   by the method, to be chosen between \code{"light"} (default) and \code{"full"}.
#'   \emph{Light} predictions produce a list including
#'   the predicted mean, standard deviation and limits of the 95\% confidence intervals at the prediction
#'   points. \emph{Full} predictions produce the same information as light ones, in addition to the
#'   training-prediction cross-covariance matrix and the prediction auto-covariance matrix.
#'
#' @return An object of class \code{"list"} containing the data structures linked to predictions. For
#'   \emph{light} predictions, the list will include the mean, standard deviation and limits of the 95\%
#'   confidence intervals at the prediction points. For \emph{full} predictions, it will include the same
#'   information, plus the training-prediction cross-covariance matrix and the prediction auto-covariance
#'   matrix.
#'
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
#'
#' @seealso \strong{*} \link[funGp]{plot.predict.fgpm} for the prediction plot of a \code{fgpm} model;
#' @seealso \strong{*} \link[funGp]{simulate,fgpm-method} for simulations based on a \code{fgpm} model;
#' @seealso \strong{*} \link[funGp]{plot.simulate.fgpm} for the simulation plot of a \code{fgpm} model.
#'
#' @examples
#' # light predictions________________________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # generating input data for prediction
#' n.pr <- 100
#' sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
#'                                 x2 = seq(0,1,length = sqrt(n.pr))))
#' fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), matrix(runif(n.pr*22), ncol = 22))
#'
#' # making predictions
#' m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
#'
#' # checking content of the list
#' summary(m1.preds)
#'
#' # ~R output:~
#' #         Length Class  Mode
#' # mean    100    -none- numeric
#' # sd      100    -none- numeric
#' # lower95 100    -none- numeric
#' # upper95 100    -none- numeric
#'
#' # plotting predictions
#' plot(m1.preds)
#'
#'
#' # comparison against true output___________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # making predictions
#' n.pr <- 100
#' sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
#'                                 x2 = seq(0,1,length = sqrt(n.pr))))
#' fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), matrix(runif(n.pr*22), ncol = 22))
#' m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
#'
#' # generating output data for validation
#' sOut.pr <- fgp_BB3(sIn.pr, fIn.pr, n.pr)
#'
#' # plotting predictions along with true output values
#' plot(m1.preds, sOut.pr)
#'
#'
#' # full predictions_________________________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # making full predictions
#' n.pr <- 100
#' sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
#'                                 x2 = seq(0,1,length = sqrt(n.pr))))
#' fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), matrix(runif(n.pr*22), ncol = 22))
#' m1.preds_f <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr, detail = "full")
#'
#' # checking content of the list
#' summary(m1.preds_f)
#'
#' # ~R output:~
#' #         Length Class  Mode
#' # mean      100  -none- numeric
#' # sd        100  -none- numeric
#' # K.tp     2500  -none- numeric
#' # K.pp    10000  -none- numeric
#' # lower95   100  -none- numeric
#' # upper95   100  -none- numeric
#'
#' # plotting predictions
#' plot(m1.preds)
#'
#' @rdname predict-methods
#' @importFrom stats qnorm
#' @importFrom methods hasArg
#' @aliases predict,fgpm-method
setMethod("predict", "fgpm",
          function(object, sIn.pr = NULL, fIn.pr = NULL, detail = c("light", "full"), ...){
              detail <- match.arg(detail)
              predict.fgpm(model = object, sIn.pr = sIn.pr, fIn.pr = fIn.pr, detail = detail)
          })

predict.fgpm <- function(model, sIn.pr, fIn.pr, detail = "light") {

  ## Added by Yves
  L <- check_new_inputs(object = model, newsIn = sIn.pr, newfIn = fIn.pr)
  ## sIn.pr <- L$newsIn
  ## fIn.pr <- L$newfIn
  ## End added by Yves

  nugget <- model@nugget

  # check validity of user inputs
  checkVal_pred_and_sim(as.list(environment()))

  if (all(model@ds > 0, model@df > 0)) { # Hybrid-input case *******************************************
    # set required data format
    sIn.pr <- as.matrix(sIn.pr)

    # project functional inputs
    f_basis <- model@f_proj@basis
    f_coefs.pr <- mapply(function(B, f) t(solve(crossprod(B), tcrossprod(t(B),f))), f_basis, fIn.pr, SIMPLIFY = FALSE)
    f_J <- lapply(f_basis, crossprod)

    # compute scalar distance matrices
    sMs.tp <- setDistMatrix_S(model@sIn, sIn.pr)
    sMs.pp <- setDistMatrix_S(sIn.pr, sIn.pr)

    # compute functional distance matrices
    fMs.tp <- setDistMatrix_F(model@f_proj@coefs, f_coefs.pr, f_J, model@kern@f_disType)
    fMs.pp <- setDistMatrix_F(f_coefs.pr, f_coefs.pr, f_J, model@kern@f_disType)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_SF(sMs.tp, sMs.pp, fMs.tp, fMs.pp,
                          model@kern@varHyp, model@kern@s_lsHyps, model@kern@f_lsHyps,
                          model@kern@kerType, model@preMats$L, model@preMats$LInvY, detail, nugget)

  } else if (model@df > 0) { # functional-input case *******************************************
    # project functional inputs
    f_basis <- model@f_proj@basis
    f_coefs.pr <- mapply(function(B, f) t(solve(crossprod(B), tcrossprod(t(B),f))), f_basis, fIn.pr, SIMPLIFY = FALSE)
    f_J <- lapply(f_basis, crossprod)

    # compute functional distance matrices
    fMs.tp <- setDistMatrix_F(model@f_proj@coefs, f_coefs.pr, f_J, model@kern@f_disType)
    fMs.pp <- setDistMatrix_F(f_coefs.pr, f_coefs.pr, f_J, model@kern@f_disType)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_F(fMs.tp, fMs.pp, model@kern@varHyp, model@kern@f_lsHyps, model@kern@kerType,
                         model@preMats$L, model@preMats$LInvY, detail, nugget)

  } else { # scalar-input case *******************************************
    # set required data format
    sIn.pr <- as.matrix(sIn.pr)

    # compute scalar distance matrices
    sMs.tp <- setDistMatrix_S(model@sIn, sIn.pr)
    sMs.pp <- setDistMatrix_S(sIn.pr, sIn.pr)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_S(sMs.tp, sMs.pp, model@kern@varHyp, model@kern@s_lsHyps, model@kern@kerType,
                         model@preMats$L, model@preMats$LInvY, detail, nugget)
  }

  # compute confidence intervals
  preds$lower95 <- preds$mean - qnorm(0.975) * preds$sd
  preds$upper95 <- preds$mean + qnorm(0.975) * preds$sd

  # _______________________________________________________________________________________________________
  # Prediction output checklist
  # _______________________________________________________________________________________________________
  # 1.  * mean ............... array (n.pr) .............. predicted mean
  # 2.  * sd ................. array (n.pr) .............. predicted standard deviation
  # 3.  * lower95 ............ array (n.pr) .............. lower bounds of 95% confidence intervals
  # 4.  * upper95 ............ array (n.pr) .............. upper bounds of 95% confidence intervals
  # 5.  * K.pp ............... matrix (n.pr x n.pr) ...... conditional covariance matrix
  # 6.  * K.tp ............... matrix (n.tr x n.pr) ...... training vs prediction cross covariance matrix
                                        # _______________________________________________________________________________________________________

    preds$mean <- drop(preds$mean)
    preds$sd <- drop(preds$sd)
    preds$lower95 <- drop(preds$lower95)
    preds$upper95 <- drop(preds$upper95)

    class(preds) <- c("predict.fgpm", "list")
    return(preds)
}
# ==========================================================================================================



# ==========================================================================================================
# Simulation based on a \code{fgpm} model
# ==========================================================================================================
#' @name simulate
#' @rdname simulate-methods
#' @importFrom stats simulate
#' @exportMethod simulate
#' @noRd
setGeneric(name = "simulate", def = function(object, nsim = 1, seed = NULL, ...) standardGeneric("simulate"))

#' @title Random sampling from a \code{fgpm} model
#' @description This method enables simulation of Gaussian process values at any given set of points
#'   based on a pre-built \code{fgpm} model. Check \code{\link{fgpm}} for information on how to create funGp models.
#'
#' @param object An object of class \linkS4class{fgpm} corresponding to the funGp model from which
#'   simulations must be performed.
#' @param nsim An optional integer indicating the number of samples to produce. Default is 1.
#' @param seed An optional value interpreted as an integer, that will be used as argument of
#'   \code{\link[base]{set.seed}} just before simulating the response values.
#' @param ... Not used.
#' @param sIn.sm An optional matrix of scalar input coordinates at which the output values should be
#'   simulated. Each column is interpreted as a scalar input variable and each row as a coordinate.
#'   Either scalar input coordinates (sIn.sm), functional input coordinates (fIn.sm), or both must be provided.
#' @param fIn.sm An optional list of functional input coordinates at which the output values should be
#'   simulated. Each element of the list is interpreted as a functional input variable. Every functional input
#'   variable should be provided as a matrix with one curve per row. Either scalar input coordinates (sIn.sm),
#'   functional input coordinates (fIn.sm), or both must be provided.
#' @param nugget.sm An optional number corresponding to a numerical nugget effect. If provided, this number
#'   is added to the main diagonal of the simulation covariance matrix in order to prevent numerical
#'   instabilities during Cholesky decomposition. A small number in the order of 1e-8 is often enough.
#'   Default is 0.
#' @param detail An optional character specifying the extent of information that should be delivered
#'   by the method, to be chosen between \code{"light"} (default)  and \code{"full"}.
#'   \emph{Light} simulations produce a matrix of
#'   simulated output values, with as many rows as requested random samples. \emph{Full} simulations produce a
#'   list with the matrix of simulated output values, along with the predicted mean, standard deviation and
#'   limits of the 95\% confidence intervals at the simulation points.
#'
#' @return An object containing the data structures linked to simulations. For \emph{light} simulations, the
#'   output will be a matrix of simulated output values, with as many rows as requested random samples.
#'   For \emph{full} simulations, the output will be a list with the matrix of simulated output values,
#'   along with the predicted mean, standard deviation and limits of the 95\% confidence intervals at the
#'   simulation points.
#'
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
#'
#' @seealso \strong{*} \link[funGp]{plot.simulate.fgpm} for the simulation plot of a \code{fgpm} model;
#' @seealso \strong{*} \link[funGp]{predict,fgpm-method} for predictions based on a \code{fgpm} model;
#' @seealso \strong{*} \link[funGp]{plot.predict.fgpm} for the prediction plot of a \code{fgpm} model.
#'
#' @examples
#' # light simulations _______________________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # generating input data for simulation
#' n.sm <- 100
#' sIn.sm <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.sm)),
#'                                 x2 = seq(0,1,length = sqrt(n.sm))))
#' fIn.sm <- list(f1 = matrix(runif(n.sm*10), ncol = 10), matrix(runif(n.sm*22), ncol = 22))
#'
#' # making light simulations
#' m1.sims_l <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm)
#'
#' # plotting light simulations
#' plot(m1.sims_l)
#'
#'
#' # full simulations ________________________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # making full simulations
#' m1.sims_f <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm, detail = "full")
#'
#' # checking content of the list
#' summary(m1.sims_f)
#'
#' # ~R output:~
#' #         Length Class  Mode
#' # sims    1000   -none- numeric
#' # mean     100   -none- numeric
#' # sd       100   -none- numeric
#' # lower95  100   -none- numeric
#' # upper95  100   -none- numeric
#'
#' # plotting full simulations in full mode
#' plot(m1.sims_f)
#'
#' # plotting full simulations in light mode
#' plot(m1.sims_f, detail = "light")
#'
#' @rdname simulate-methods
#' @aliases simulate,fgpm-method
setMethod("simulate", "fgpm",
          function(object, nsim = 1, seed = NULL, sIn.sm = NULL, fIn.sm = NULL,
                   nugget.sm = 0, detail = c("light", "full"), ...) {
              detail <- match.arg(detail)
            simulate.fgpm(model = object, nsim = nsim, seed = seed, sIn.sm = sIn.sm, fIn.sm = fIn.sm,
                           nugget.sm = nugget.sm, detail = detail)
          })

simulate.fgpm <- function(model, nsim, seed, sIn.sm, fIn.sm, nugget.sm = 10^-8, detail) {
  # check validity of user inputs
  checkVal_pred_and_sim(as.list(environment()))

  ## Added by Yves
  L <- check_new_inputs(object = model, newsIn = sIn.sm, newfIn = fIn.sm)
  ## sIn.sm <- L$newsIn
  ## fIn.sm <- L$newfIn
  ## End added by Yves


  # check which type of model it is
  if (all(model@ds > 0, model@df > 0)) { # Hybrid-input case *******************************************
    # set required data format
    sIn.sm <- as.matrix(sIn.sm)

    # project functional inputs
    f_basis <- model@f_proj@basis
    f_coefs.sm <- mapply(function(B, f) t(solve(crossprod(B), tcrossprod(t(B),f))), f_basis, fIn.sm, SIMPLIFY = FALSE)
    f_J <- lapply(f_basis, crossprod)

    # compute scalar distance matrices
    sMs.ts <- setDistMatrix_S(model@sIn, sIn.sm)
    sMs.ss <- setDistMatrix_S(sIn.sm, sIn.sm)

    # compute functional distance matrices
    fMs.ts <- setDistMatrix_F(model@f_proj@coefs, f_coefs.sm, f_J, model@kern@f_disType)
    fMs.ss <- setDistMatrix_F(f_coefs.sm, f_coefs.sm, f_J, model@kern@f_disType)

    # make simulations based on the Gaussian Conditioning Theorem
    sims <- makeSims_SF(sMs.ts, sMs.ss, fMs.ts, fMs.ss,
                        model@kern@varHyp, model@kern@s_lsHyps, model@kern@f_lsHyps,
                        model@kern@kerType, model@preMats$L, model@preMats$LInvY, nsim, nugget.sm, detail, seed)

  } else if (model@df > 0) { # functional-input case *******************************************
    # project functional inputs
    f_basis <- model@f_proj@basis
    f_coefs.sm <- mapply(function(B, f) t(solve(crossprod(B), tcrossprod(t(B),f))), f_basis, fIn.sm, SIMPLIFY = FALSE)
    f_J <- lapply(f_basis, crossprod)

    # compute functional distance matrices
    fMs.ts <- setDistMatrix_F(model@f_proj@coefs, f_coefs.sm, f_J, model@kern@f_disType)
    fMs.ss <- setDistMatrix_F(f_coefs.sm, f_coefs.sm, f_J, model@kern@f_disType)

    # make simulations based on the Gaussian Conditioning Theorem
    sims <- makeSims_F(fMs.ts, fMs.ss, model@kern@varHyp, model@kern@f_lsHyps, model@kern@kerType,
                       model@preMats$L, model@preMats$LInvY, nsim, nugget.sm, detail, seed)

  } else { # scalar-input case *******************************************
    # set required data format
    sIn.sm <- as.matrix(sIn.sm)

    # compute scalar distance matrices
    sMs.ts <- setDistMatrix_S(model@sIn, sIn.sm)
    sMs.ss <- setDistMatrix_S(sIn.sm, sIn.sm)

    # make simulations based on the Gaussian Conditioning Theorem
    sims <- makeSims_S(sMs.ts, sMs.ss, model@kern@varHyp, model@kern@s_lsHyps, model@kern@kerType,
                       model@preMats$L, model@preMats$LInvY, nsim, nugget.sm, detail, seed)
  }

  # if detail == 'full', confidence intervals at simulation points are provided,
  # else the sims list is dropped to a matrix with the observations only
    if (detail == "full") {
        sims$mean <- drop(sims$mean)
        sims$sd <- drop(sims$sd)

    # compute confidence intervals
    sims$lower95 <- sims$mean - qnorm(0.975) * sims$sd
    sims$upper95 <- sims$mean + qnorm(0.975) * sims$sd
  } else {
    sims <- sims$sims
  }

  # _______________________________________________________________________________________________________
  # Simulation output checklist
  # _______________________________________________________________________________________________________
  # 1.  * sims ............... matrix (nsim x n.sm) ...... simulated output
  # 2.  * mean ............... array (n.sm) .............. predicted mean
  # 3.  * sd ................. array (n.sm) .............. predicted standard deviation
  # 4.  * lower95 ............ array (n.sm) .............. lower bounds of 95% confidence intervals
  # 5.  * upper95 ............ array (n.sm) .............. upper bounds of 95% confidence intervals
  # _______________________________________________________________________________________________________


    class(sims) <- c("simulate.fgpm", "list")
    return(sims)


  return(sims)
}
# ==========================================================================================================



# ==========================================================================================================
# Updating of a \code{fgpm} model
# ==========================================================================================================
#' @name update
#' @rdname update-methods
#' @importFrom stats update
#' @exportMethod update
#' @noRd
setGeneric(name = "update", def = function(object, ...) standardGeneric("update"))

#' @title Easy update of \code{fgpm} models
#' @description This method enables the update of data or hyperparameters of a \code{fgpm} model.
#'   It corresponds to an object of the class \linkS4class{fgpm}. The method allows addition, subtraction
#'   and substitution of data points, as well as substitution and re-estimation of hyperparameters.
#'
#' @param object An object of class \linkS4class{fgpm} corresponding to the funGp model to update.
#' @param ... Not used.
#' @param sIn.nw An optional matrix of scalar input values to be added to the model. Each column must match
#'   an input variable and each row a scalar coordinate.
#' @param fIn.nw An optional list of functional input values to be added to the model. Each element of the
#'   list must be a matrix containing the set of curves corresponding to one functional input.
#' @param sOut.nw An optional vector (or 1-column matrix) containing the values of the scalar output at the
#'   new input points.
#' @param sIn.sb An optional matrix of scalar input values to be used as substitutes of other scalar input
#'   values already stored in the model. Each column must match an input variable and each row a coordinate.
#' @param fIn.sb An optional list of functional input values to be added to the model. Each element of the
#'   list must be a matrix containing the set of curves corresponding to one functional input.
#' @param sOut.sb An optional vector (or 1-column matrix) containing the values of the scalar output at the
#'   substituting input points.
#' @param ind.sb An optional numeric array indicating the indices of the input and output points stored in
#'   the model, that should be replaced by the values specified through sIn.sb, fIn.sb and/or sOut.sb.
#' @param ind.dl An optional numeric array indicating the indices of the input and output points stored in
#'   the model that should be deleted.
#' @param var.sb An optional number indicating the value that should be used to substitute the current
#'   variance parameter of the model.
#' @param ls_s.sb An optional numerical array indicating the values that should be used to substitute the
#'   current length-scale parameters for the scalar inputs of the model.
#' @param ls_f.sb An optional numerical array indicating the values that should be used to substitute the
#'   current length-scale parameters for the functional inputs of the model.
#' @param var.re An optional boolean indicating whether the variance parameter should be re-estimated.
#'   Default is FALSE.
#' @param ls_s.re An optional boolean indicating whether the length-scale parameters of the scalar inputs
#'   should be re-estimated. Default is FALSE.
#' @param ls_f.re An optional boolean indicating whether the length-scale parameters of the functional
#'   inputs should be re-estimated. Default is FALSE.
#' @param trace An optional boolean indicating whether a summary update should be displayed. Default is TRUE.
#'
#' @return An object of class \linkS4class{fgpm} representing the updated funGp model.
#'
#' @details
#' The arguments listed above enable the completion of the following updating tasks:
#' \itemize{
#'  \item \strong{Deletion} of data points: ind.dl;
#'  \item \strong{Addition} of data points: sIn.nw, fIn.nw, sOut.nw;
#'  \item \strong{Substitution} of data points: sIn.sb, fIn.sb, sOut.sb, ind.sb;
#'  \item \strong{Substitution} of hyperparameters: var.sb, ls_s.sb, ls_f.sb;
#'  \item \strong{Re-estimation} of hyperparameters: var.re, ls_s.re, ls_f.re.
#' }
#'
#' All the arguments listed above are optional since any of these tasks can be requested without need to
#' request any of the other tasks. In fact, most of the arguments can be used even if the other
#' arguments related to the same task are not. For instance, the re-estimation of the variance can be
#' requested via var.re without requiring re-estimation of the scalar or functional length-scale
#' parameters. The only two exceptions are: (i) for data addition, the new output sOut.nw should always
#' be provided and the new input points should correspond to the set of variables already stored in the
#' \linkS4class{fgpm} object passed for update; and (ii) for data substitution, the argument ind.sb is
#' always mandatory.
#'
#' @details
#' \strong{Conflicting task combinations:}
#' \itemize{
#'  \item Data points deletion and substitution;
#'  \item Substitution and re-estimation of the same hyperparameter.
#' }
#'
#' @details
#' Note that the parameters of the model will not be updated after modifying the model unless explicitly
#' requested through the var.re, ls_s.re and ls_f.re arguments. If, for instance, some points are added
#' to the model without requesting parameter re-estimation, the new data will be included in the
#' training-training and training-prediction covariance matrices, but the hyperparameters will not
#' be updated. This allows to make updates in the data that might help to improve predictions,
#' without the immediate need to perform a training procedure that could be time consuming. At any later
#' time, the user is allowed to request the re-estimation of the hyperparameters, which will make
#' the model fully up to date.
#'
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
#'
#' @seealso \strong{*} \link[funGp]{fgpm} for creation of a funGp model;
#' @seealso \strong{*} \link[funGp]{predict,fgpm-method} for predictions based on a \code{fgpm} model;
#' @seealso \strong{*} \link[funGp]{simulate,fgpm-method} for simulations based on a \code{fgpm} model.
#'
#' @examples
#' # deletion and addition of data points_____________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # deleting two points
#' ind.dl <- sample(1:m1@n.tot, 2)
#' m1up <- update(m1, ind.dl = ind.dl)
#'
#' # adding five points
#' n.nw <- 5
#' sIn.nw <- matrix(runif(n.nw * m1@ds), nrow = n.nw)
#' fIn.nw <- list(f1 = matrix(runif(n.nw*10), ncol = 10), f2 = matrix(runif(n.nw*22), ncol = 22))
#' sOut.nw <- fgp_BB3(sIn.nw, fIn.nw, n.nw)
#' m1up <- update(m1, sIn.nw = sIn.nw, fIn.nw = fIn.nw, sOut.nw = sOut.nw)
#'
#'
#' # substitution of data points______________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # generating substituting input data for updating
#' n.sb <- 2
#' sIn.sb <- matrix(runif(n.sb * m1@ds), nrow = n.sb)
#' fIn.sb <- list(f1 = matrix(runif(n.sb*10), ncol = 10), f2 = matrix(runif(n.sb*22), ncol = 22))
#'
#' # generating substituting output data for updating
#' sOut.sb <- fgp_BB3(sIn.sb, fIn.sb, n.sb)
#'
#' # generating indices for substitution
#' ind.sb <- sample(1:(m1@n.tot), n.sb)
#'
#' # updating all, the scalar inputs, functional inputs and the outputs
#' m1up <- update(m1, sIn.sb = sIn.sb, fIn.sb = fIn.sb, sOut.sb = sOut.sb, ind.sb = ind.sb)
#'
#' # updating only some of the data structures
#' m1up1 <- update(m1, sIn.sb = sIn.sb, ind.sb = ind.sb) # only the scalar inputs
#' m1up2 <- update(m1, sOut.sb = sOut.sb, ind.sb = ind.sb) # only the outputs
#' m1up3 <- update(m1, sIn.sb = sIn.sb, sOut.sb = sOut.sb, ind.sb = ind.sb) # the scalar inputs
#'                                                                          # and the outputs
#'
#'
#' # substitution of hyperparameters__________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # defining hyperparameters for substitution
#' var.sb <- 3
#' ls_s.sb <- c(2.44, 1.15)
#' ls_f.sb <- c(5.83, 4.12)
#'
#' # updating the model
#' m1up <- update(m1, var.sb = var.sb, ls_s.sb = ls_s.sb, ls_f.sb = ls_f.sb)
#'
#'
#' # re-estimation of hyperparameters_________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # re-estimating the hyperparameters
#' m1up <- update(m1, var.re = TRUE) # only the variance
#' m1up <- update(m1, ls_s.re = TRUE) # only the scalar length-scale parameters
#' m1up <- update(m1, ls_s.re = TRUE, ls_f.re = TRUE) # all length-scale parameters
#' m1up <- update(m1, var.re = TRUE, ls_s.re = TRUE, ls_f.re = TRUE) # all hyperparameters
#'
#' @rdname update-methods
#' @importFrom utils tail
#' @aliases update,fgpm-method
setMethod("update", "fgpm",
          function(object, sIn.nw = NULL, fIn.nw = NULL, sOut.nw = NULL,
                   sIn.sb = NULL, fIn.sb = NULL, sOut.sb = NULL, ind.sb = NULL,
                   ind.dl = NULL, var.sb = NULL, ls_s.sb = NULL, ls_f.sb = NULL,
                   var.re = FALSE, ls_s.re = FALSE, ls_f.re = FALSE, trace = TRUE, ...) {
            update.fgpm(model = object, sIn.nw = sIn.nw, fIn.nw = fIn.nw, sOut.nw = sOut.nw,
                        sIn.sb = sIn.sb, fIn.sb = fIn.sb, sOut.sb = sOut.sb, ind.sb = ind.sb,
                        ind.dl = ind.dl,
                        var.sb = var.sb, ls_s.sb = ls_s.sb, ls_f.sb = ls_f.sb,
                        var.re = var.re, ls_s.re = ls_s.re, ls_f.re = ls_f.re,
                        trace = trace)

          })

update.fgpm <- function(model, sIn.nw, fIn.nw, sOut.nw, sIn.sb, fIn.sb, sOut.sb, ind.sb, ind.dl,
                         var.sb, ls_s.sb, ls_f.sb, var.re, ls_s.re, ls_f.re, trace) {
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
                 "var substitution", "scalar length-scale substitution",
                 "functional length-scale substitution",
                 "var re-estimation", "scalar length-scale re-estimation",
                 "functional length-scale re-estimation")

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
    modelup@howCalled <- model@howCalled
    modelup@n.tr <- model@n.tr
    modelup@convergence <- model@convergence
    modelup@NegLogLike <- model@NegLogLike
    cptasks <- c(cptasks, 1)
  }
  if (subInOut & !(2 %in% dptasks)) {
    modelup <- upd_subData(model = modelup, sIn.sb = sIn.sb, fIn.sb = fIn.sb,
                           sOut.sb = tryCatch(as.matrix(sOut.sb), error = function(e) sOut.sb), ind.sb = ind.sb,
                           remake = all(!newInOut, !subHypers, !reeHypers))
    modelup@howCalled <- model@howCalled
    modelup@convergence <- model@convergence
    modelup@NegLogLike <- model@NegLogLike
    cptasks <- c(cptasks, 2)
  }
  if (newInOut) {
    modelup <- upd_add(model = modelup, sIn.nw = sIn.nw, fIn.nw = fIn.nw, sOut.nw = as.matrix(sOut.nw),
                       remake = all(!subHypers, !reeHypers))
    modelup@howCalled <- model@howCalled
    modelup@n.tr <- model@n.tr
    modelup@convergence <- model@convergence
    modelup@NegLogLike <- model@NegLogLike
    cptasks <- c(cptasks, 3)
  }
  if (subHypers & any(!(c(4,5,6) %in% dptasks))) {
    modelup <- upd_subHypers(model = modelup, var.sb = var.sb, ls_s.sb = ls_s.sb, ls_f.sb = ls_f.sb)
    modelup@howCalled <- model@howCalled
    modelup@n.tr <- model@n.tr
    modelup@convergence <- model@convergence
    modelup@NegLogLike <- model@NegLogLike
    if (!is.null(var.sb) & !(4 %in% dptasks)) cptasks <- c(cptasks, 4)
    if (!is.null(ls_s.sb) & !(5 %in% dptasks)) cptasks <- c(cptasks, 5)
    if (!is.null(ls_f.sb) & !(6 %in% dptasks)) cptasks <- c(cptasks, 6)
  }
  if (reeHypers & any(!(c(7,8,9) %in% dptasks))) {
    modelup <- upd_reeHypers(model = modelup, var.re = var.re, ls_s.re = ls_s.re, ls_f.re = ls_f.re)
    modelup@howCalled <- model@howCalled
    if (isTRUE(var.re) & !(7 %in% dptasks)) cptasks <- c(cptasks, 7)
    if (isTRUE(ls_s.re) & !(8 %in% dptasks)) cptasks <- c(cptasks, 8)
    if (isTRUE(ls_f.re) & !(9 %in% dptasks)) cptasks <- c(cptasks, 9)
  }
  # ----------------------------------------------------

  # print update summary (only if trace is enabled)
  if (trace) {
    # ----------------------------------------------------
    if (length(cptasks) > 0) { # list of complete tasks if there is any
      cat("--------------\n")
      cat("Update summary\n")
      cat("--------------\n")

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
  }
  # ----------------------------------------------------

  return(modelup)
}

## ==============================================================================
## summary method. Simple copy of 'show', at least for now.
## ==============================================================================
##' @description Display the structure of a \code{fgpm}
##'     object and the value of the parameters (variance and length-scales).
##'
##' @title Summary method for \code{fgpm} objects
##' @param object An \code{fgpm} object.
##' @param ... Not used yet.
##' @method summary fgpm
##'
##' @note This method is actually identical to the \code{show} method
##'     for this class which is called when the name of the object is
##'     entered in an interactive session.
##'
##' @examples
##' m <- xm@model
##' class(m)
##' summary(m)
##' m
setMethod("summary", "fgpm", function(object, ...) show(object, ...))

