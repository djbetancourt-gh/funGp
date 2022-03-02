# ==========================================================================================================
# S4 class for funGp model selection data structures
# ==========================================================================================================
#' @title S4 class for funGp model selection data structures
#' @description This is the formal representation of the assembly of data structures delivered by the model
#' selection routines in the \link[=funGp-package]{funGp package}. Gaussian process models are useful
#' statistical tools in the modeling of complex input-output relationships. An Xfgpm object contains the
#' trace of an optimization process, conducted to build Gaussian process models of outstanding performance.
#' \itemize{
#'  \item{\strong{Main methods}}{\cr
#'        \link[funGp]{fgpm_factory}: structural optimization of funGp models \cr
#'  }
#'  \item{\strong{Plotters}}{\cr
#'        \link[funGp]{plot,Xfgpm-method}: plot of the evolution of the algorithm with \code{which = "evolution"}
#'        or of the absolute and relative quality of the optimized model with \code{which = "diag"}
#'  }
#' }
#'
#' @slot factoryCall Object of class \code{"\linkS4class{factoryCall}"}. User call reminder.
#' @slot model Object of class \code{"\linkS4class{fgpm}"}. Model selected by the heuristic structural
#'   optimization algorithm.
#' @slot stat Object of class \code{"character"}. Performance measure optimized to select the model. To be
#'   set from "Q2loocv", "Q2hout".
#' @slot fitness Object of class \code{"numeric"}. Value of the performance measure for the selected model.
#' @slot structure Object of class \code{"data.frame"}. Structural configuration of the selected model.
#' @slot log.success Object of class \code{"\linkS4class{antsLog}"}. Record of models successfully
#'   evaluated during the structural optimization. It contains the structural configuration both in
#'   data.frame and \code{"\linkS4class{modelCall}"} format, along with the fitness of each model. The
#'   models are sorted by fitness, starting with the best model in the first position.
#' @slot log.crashes Object of class \code{"\linkS4class{antsLog}"}. Record of models crashed during the
#'   structural optimization. It contains the structural configuration of each model, both in data.frame
#'   and \code{"\linkS4class{modelCall}"} format.
#' @slot n.solspace Object of class \code{"numeric"}. Number of possible structural configurations for
#'   the optimization instance resolved.
#' @slot n.explored Object of class \code{"numeric"}. Number of structural configurations successfully
#'   evaluated by the algorithm.
#' @slot details Object of class \code{"list"}. Further information about the parameters of the ant colony
#'   optimization algorithm and the evolution of the fitness along the iterations.
#' @slot sIn An object of class \code{"matrix"} containing a copy of
#'     the provided scalar inputs.
#' @slot fIn An object of class \code{"list"} containing a copy of
#'     the provided functional inputs.
#' @slot sOut An object of class \code{"matrix"} containing a copy of the provided output.
#' @section Useful material:
#' \itemize{
#'  \item{\strong{Manual}}{
#'  \href{https://hal.archives-ouvertes.fr/hal-02536624}{
#'  Gaussian Process Regression for Scalar and Functional Inputs with funGp - The in-depth tour
#'  }}
#' }
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @include 1_fgpm_Class.R
#' @include 3_ant_admin.R
#' @include 8_outilsCode.R
#'
#' @rdname Xfgpm-class
#' @export
setClass("Xfgpm",
         representation(
           factoryCall = "factoryCall",    # distance type. To be chosen from {"scalar", "functional"}
           model = "fgpm",                 # kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}
           stat = "character",             # search method
           fitness = "numeric",            # model fitness
           structure = "data.frame",       # model fitness
           log.success = "antsLog",        # search method
           log.crashes = "antsLog",        # search method
           n.solspace = "numeric",         # search method
           n.explored = "numeric",         # search method
           details = "list",               # search method
           sIn = "matrix",
           fIn = "list",
           sOut = "matrix"
         ),
         validity = function(object) {TRUE})



# ==========================================================================================================
# Printing of an XfunGp model
# ==========================================================================================================
#' @rdname show-methods
#' @aliases show,Xfgpm-method
setMethod("show", "Xfgpm", function(object) show.Xfgpm(object))

show.Xfgpm <- function(object) {
  cat("Structural optimization_________________\n")

  cat(paste("stat:", object@stat, "\n"))
  cat(paste("value:", format(object@fitness, digits = 3, nsmall = 4), "\n"))
  cat(paste("n.explored:", object@n.explored, "\n"))

  cat("For selected structure: object@structure\n")
  cat("For log of success: object@log.success\n")
  cat("________________________________________\n")
}
# ==========================================================================================================



# ==========================================================================================================
# Structural optimization of a funGp model
# ==========================================================================================================
#' @title Structural optimization of Gaussian process models
#' @description This function enables the smart exploration of the solution space of potential structural
#'   configurations of a funGp model, and the consequent selection of a high quality configuration. funGp
#'   currently relies on an ant colony based algorithm to perform this task. The algorithm defines the
#'   solution space based on the levels of each structural parameter currently available in the
#'   \link[funGp]{fgpm} function, and performs as smart exploration of it. More details on the algorithm are
#'   provided in a dedicated
#'   \href{https://hal.archives-ouvertes.fr/hal-02532713}{technical report}.
#'   funGp might evolve in the future to include improvements in the current algorithm or alternative
#'   solution methods.
#'
#' @param sIn An optional matrix of scalar input values to train the model. Each column must match an input
#'   variable and each row a training point. Either scalar input coordinates (sIn), functional input
#'   coordinates (fIn), or both must be provided.
#' @param fIn An optional list of functional input values to train the model. Each element of the list must
#'   be a matrix containing the set of curves corresponding to one functional input. Either scalar input
#'   coordinates (sIn), functional input coordinates (fIn), or both must be provided.
#' @param sOut A vector (or 1-column matrix) containing the values of the scalar output at the specified
#'   input points.
#' @param ind.vl An optional numerical matrix specifying which points in the three structures above should be
#'   used for training and which for validation. If provided, the optimization will be conducted in terms of
#'   the hold-out Q2, which comes from training the model with a subset of the points, and then estimate the
#'   prediction error in the remaining points. In that case, each column of \emph{ind.vl} will be interpreted
#'   as one validation set, and the multiple columns will imply replicates. In the simplest case,
#'   \emph{ind.vl} will be a one-column matrix or simply an array, meaning that a simple replicate should be
#'   used for each model configuration explored. If not provided, the optimization will be conducted in terms
#'   of the leave-one-out cross-validation Q2, which for a total number of n observations, comes from training
#'   the model n times, each using n-1 points for training and the remaining one for validation. This procedure
#'   is typically costly due to the large number of hyperparameter optimizations that should be conducted,
#'   nonetheless, fgpm_factory implements the virtual equations introduced by Dubrule (1983) for Gaussian
#'   processes, which require a single hyperparameter optimization. See the reference below for more details.
#' @param ctraints An optional list specifying the constraints of the structural optimization problem. Valid
#'   entries for this list are: \cr\cr
#'   \strong{*}\emph{s_keepOn}: a numerical array indicating the scalar inputs that should remain active in the
#'     model. It should contain the index of the columns of sIn corresponding to the inputs to keep active. \cr\cr
#'   \strong{*}\emph{f_keepOn}: a numerical array indicating the functional inputs that should remain active in
#'     the model. It should contain the index of the elements of fIn corresponding to the inputs to keep active. \cr\cr
#'   \strong{*}\emph{f_disTypes}: a list specifying the set of distances that should be tested for some
#'     functional inputs. The values should be taken from the possibilities offered by the \link[funGp]{fgpm}
#'     function for the argument \emph{f_disType} therein. Valid choices at this time are "L2_bygroup" and
#'     "L2_byindex". Each element of the list should receive as name the index of a functional input variable,
#'     and should contain an array of strings with the name of the distances allowed for this input. All the
#'     available distances will be tried for any functional input not included in the list. \cr\cr
#'   \strong{*}\emph{f_fixDims}: a two-row matrix specifying a particular projection dimension for some
#'     functional inputs. For each input, the value should be a number between 0 and its original dimension,
#'     with 0 denoting no projection. The first row of the matrix should contain the index of each input, and
#'     the second row should contain the corresponding dimensions. All the possible dimensions will be tried
#'     for any functional input not included in the matrix (unless affected by the \emph{f_maxDims} argument
#'     below). \cr\cr
#'   \strong{*}\emph{f_maxDims}: a two-row matrix specifying the largest projection dimension for some
#'     functional inputs. For each input, the value should be a number between 1 and its original dimension.
#'     The first row of the matrix should contain the index of each input, and the second row should contain
#'     the corresponding largest dimensions. All the possible dimensions will be tried for any functional input
#'     not included in the matrix (unless affected by the \emph{f_fixDims} argument above). \cr\cr
#'   \strong{*}\emph{f_basTypes}: a list specifying the set of basis families that should be tested for some
#'     functional inputs. The values should be taken from the possibilities offered by the \link[funGp]{fgpm}
#'     function for the argument \emph{f_basType} therein. Valid choices at this time are "B-splines" and "PCA".
#'     Each element of the list should receive as name the index of a functional input variable, and should
#'     contain an array of strings with the name of the distances allowed for this input. All the available
#'     basis families will be tried for any functional input not included in the list. \cr\cr
#'   \strong{*}\emph{kerTypes}: an array of strings specifying the kernel functions allowed to be tested. The
#'     values should be taken from the possibilities offered by the \link[funGp]{fgpm} function for the argument
#'     \emph{kerType} therein. Valid choices at this time are "gauss", "matern5_2" and "matern3_2". If not
#'     provided, all the available kernel functions will be tried.
#' @param setup An optional list indicating the value for some parameters of the structural optimization
#'   algorithm. The ant colony optimization algorithm available at this time allows the following entries: \cr\cr
#'   \strong{Initial pheromone load}\cr\cr
#'     \strong{*}\emph{tao0}: a number indicating the initial pheromone load on links pointing out to the
#'       selection of a distance type, a projection basis or a kernel type. Default is 0.1. \cr\cr
#'     \strong{*}\emph{dop.s}: a number controlling how likely it is to activate a scalar input. It operates on a
#'       relation of the type \eqn{A  = dop.s * I}, where \emph{A} is the initial pheromone load of links
#'       pointing out to the activation of scalar inputs and \emph{I} is the initial pheromone load of links
#'       pointing out to their inactivation. Default is 1. \cr\cr
#'     \strong{*}\emph{dop.f}: analogous to \emph{dop.s} for functional inputs. Default is 1. \cr\cr
#'     \strong{*}\emph{delta.f and dispr.f}: two numbers used as shape parameters for the regularization
#'       function that determines the initial pheromone values on the links connecting the L2_byindex distance
#'       with the projection dimension. Default are 2 and 1.4, respectively. \cr\cr
#'   \strong{Local pheromone update}\cr\cr
#'     \strong{*}\emph{rho.l}: a number specifying the pheromone evaporation rate. Default is 0.1. \cr\cr
#'   \strong{Global pheromone update}\cr\cr
#'     \strong{*}\emph{u.gbest}: a boolean indicating if at each iteration, the pheromone load on the links
#'       of the best ant of the whole trial should be reinforced. Default is FALSE. \cr\cr
#'     \strong{*}\emph{n.ibest}: a number indicating how many top ants of each iteration should be used for
#'       pheromone reinforcement. Default is 1. \cr\cr
#'     \strong{*}\emph{rho.g}: a number specifying the learning reinforcement rate. Default is 0.1. \cr\cr
#'   \strong{Population factors}\cr\cr
#'     \strong{*}\emph{n.iter}: a number specifying the amount of iterations of the algorithm. Default is 15. \cr\cr
#'     \strong{*}\emph{n.pop}: a number specifying the amount of ants per iteration; each ant corresponds to one
#'       structural configuration for the model. Default is 10. \cr\cr
#'   \strong{Bias strength}\cr\cr
#'     \strong{*}\emph{q0}: ants use one of two rules to select their next node at each step. The first rule leads
#'       the ant through the link with higher pheromone load; the second rule works based on probabilities which
#'       are proportional to the pheromone load on the feasible links. The ants will randomly chose one of the two
#'       rules at each time. They will opt for rule 1 with probability \emph{q0}. Default is 0.95.
#' @param time.lim An optional number specifying a time limit in seconds to be used as stopping condition for the
#'   structural optimization.
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
#' @param par.clust An optional parallel processing cluster created with the \code{\link[parallel]{makeCluster}}
#'   function of the \link[=parallel]{parallel package}. If not provided, structural configurations are evaluated in
#'   sequence.
#' @param pbars An optional boolean indicating if progress bars should be displayed.
#'
#' @return An object of class \linkS4class{Xfgpm} containing the data structures linked to the structural optimization
#'   of a funGp model. It includes as the main component, an object of class \linkS4class{fgpm} corresponding to the
#'   optimized model. It is accessible through the \code{@@model} slot of the Xfgpm object.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
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
#' @references Dubrule, O. (1983),
#' "Cross validation of kriging in a unique neighborhood".
#' \emph{Journal of the International Association for Mathematical Geology}, \strong{15},  687-699.
#' \href{https://link.springer.com/article/10.1007/BF01033232}{[MG]}
#'
#'
#'
#' @seealso \strong{*} \link[funGp]{plot,Xfgpm-method} for a call to \link[funGp]{plotEvol} with \code{which = "evolution"}
#'        or to \link[funGp]{plotX} with \code{which = "diag"};
#' @seealso \strong{*} \link[funGp]{plotEvol} for a plot of the evolution of the model selection algorithm in fgpm_factory;
#' @seealso \strong{*} \link[funGp]{plotX} for diagnostic plots for a fgpm_factory output and selected model;
#' @seealso \strong{*} \link[funGp]{get_active_in} for post-processing of input data structures following a fgpm_factory call;
#' @seealso \strong{*} \link[funGp]{predict} for predictions based on a funGp model;
#' @seealso \strong{*} \link[funGp]{simulate} for simulations based on a funGp model;
#' @seealso \strong{*} \link[funGp]{update} for post-creation updates on a funGp model.
#'
#' @examples
#' # Data with precalculated Xfgpm objects are already available
#' # Optimized model structure for fgp_BB7 black-box function with standard parameters
#' # assessing the quality of the model
#' # in the absolute and also w.r.t. the other explored models
#' plot(xm, which="diag")
#'
#' # checking the evolution of the algorithm
#' plot(xm, which="evol")
#'
#' # Summary of the tested configurations
#' summary(xm)
#'
#' # checking the log of crashed iterations
#' print(xm@log.crashes)
#'
#' # building the model with the default fgpm arguments to compare
#' set.seed(100)
#' n.tr <- 32
#' x1 <- x2 <- x3 <- x4 <- x5 <- seq(0,1,length = n.tr^(1/5))
#' sIn <- expand.grid(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
#' fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
#' f2 <- matrix(runif(n.tr * 22), ncol = 22))
#' sOut <- fgp_BB7(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#' plot(m1) # plotting the model
#'
#' # improving performance with more iterations_______________________________________________
#' plot(xm25, which="evol")
#' plot(xm25, which="diag")
#'
#' # custom solution space____________________________________________________________________
#' plot(xmc, which="evol")
#' plot(xmc, which="diag")
#'
#' # verifying constraints with the log of some successfully built models
#' summary(xmc)
#'
#' # custom heuristic parameters______________________________________________________________
#' # verifying heuristic setup through the details of the Xfgpm object
#' unlist(xmh@details$param)
#'
#' # stopping condition based on time_________________________________________________________
#' summary(xms)
#'
#' \dontrun{
#' # parallelization in the model factory_____________________________________________________
#' # generating input and output data
#' set.seed(100)
#' n.tr <- 243
#' sIn <- expand.grid(x1 = seq(0,1,length = n.tr^(1/5)), x2 = seq(0,1,length = n.tr^(1/5)),
#'                    x3 = seq(0,1,length = n.tr^(1/5)), x4 = seq(0,1,length = n.tr^(1/5)),
#'                    x5 = seq(0,1,length = n.tr^(1/5)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB7(sIn, fIn, n.tr)
#'
#' # calling fgpm_factory in parallel
#' cl <- parallel::makeCluster(2)
#' xm.par <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut, par.clust = cl) #  (~260 seconds)
#' parallel::stopCluster(cl)
#'
#' # NOTE: in order to provide progress bars for the monitoring of time consuming processes
#' #       ran in parallel, funGp relies on the doFuture and future packages. Parallel processes
#' #       suddenly interrupted by the user tend to leave corrupt connections. This problem is
#' #       originated outside funGp, which limits our control over it. On section 4.1 of the manual
#' #       of funGp, we provide a temporary solution to the issue and we remain attentive in
#' #       case it appears a more elegant way to handle it or a manner to suppress it.
#' #
#' #       funGp manual: https://hal.archives-ouvertes.fr/hal-02536624
#' }
#'
#' @importFrom methods new
#' @importFrom microbenchmark microbenchmark
#' @export
fgpm_factory <- function(sIn = NULL, fIn = NULL, sOut = NULL, ind.vl = NULL,
                          ctraints = list(), setup = list(), time.lim = Inf,
                          nugget = 1e-8, n.starts = 1, n.presample = 20,
                          par.clust = NULL, pbars = interactive()) {

  # launch timer
  time.str <- Sys.time()

  # extra arguments for model call
  extargs <- list(nugget = nugget, n.starts = n.starts, n.presample = n.presample)

  # define solution space based on user inputs
  solspace <- setSpace(sIn, fIn, ctraints)

  if (!is.null(ind.vl)) {
    ind.vl <- as.matrix(ind.vl)
    stat <- paste("Q2hout (", (nrow(sOut) - nrow(ind.vl)), ".", nrow(ind.vl), ".", ncol(ind.vl), ")", sep = "")
  } else {
    stat <- "Q2loocv"
  }

  # optimize model structure
  opt <- master_ACO(sIn, fIn, sOut, ind.vl, solspace, setup, extargs, time.str, time.lim, pbars, par.clust)
  X.model <- new("Xfgpm")
  X.model@factoryCall@string <- gsub("^ *|(?<= ) | *$", "", paste0(deparse(match.call()), collapse = " "), perl = TRUE)
  X.model@model <- opt$model
  X.model@stat <- stat
  X.model@fitness <- opt$b.fitness
  X.model@structure <- opt$sol.vec
  X.model@log.success <- opt$log.suc
  X.model@log.crashes <- opt$log.cra
  X.model@n.solspace <- getSpacesize(solspace$sp.user)
  X.model@n.explored <- nrow(opt$log.suc@sols)
    X.model@details <- opt$all.details
    ## XXXY warn if these are not matrices, but data frames?
    X.model@sIn <- as.matrix(sIn)
    X.model@fIn <- fIn
    X.model@sOut <- as.matrix(sOut)
    return(X.model)
}
# ==========================================================================================================



# ==========================================================================================================
# ==========================================================================================================
setSpace <- function(sIn, fIn, ctraints) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  # recover individual ctraints from user inputs
  s_keepOn <- ctraints$s_keepOn
  f_keepOn <- ctraints$f_keepOn
  f_fixDims <- ctraints$f_fixDims
  f_maxDims <- ctraints$f_maxDims
  f_disTypes <- ctraints$f_disTypes
  f_basTypes <- ctraints$f_basTypes
  kerTypes <- ctraints$kerTypes

  # define template with base structures and values for a full experiment
  # s.state: 0 = free, 1 = fixed
  # f.state: 0 = free, 1 = fixed
  # f.dims: seq. of potential dimensions f/e variable. 0 = no projection
  # f.dist: all available distances at the time
  # f.bas: all available basis families for functions at the time
  # k.type: all available kernels at the time
  s.state <- rep(0, ds)
  f.state <- rep(0, df)
  f.dims <- lapply(sapply(fIn, ncol), function(k) 0:k)
  f.dist <- rep(list(c("L2_bygroup", "L2_byindex")), df)
  f.bas <- rep(list(c("B-splines", "PCA")), df)
  k.type <- c("gauss", "matern5_2", "matern3_2")
  sp.base <- list(s.state = s.state, f.state = f.state, f.dims = f.dims, f.dist = f.dist, f.bas = f.bas, k.type = k.type)

  # modify the solution space if the user has specified any constraint
  if (length(ctraints) > 0) {
    # update state of scalar inputs
    if (!is.null(s_keepOn)) s.state[s_keepOn] <- 1

    # update state of functional inputs
    if (!is.null(f_keepOn)) f.state[f_keepOn] <- 1

    # update the set of potential dimensions for functional inputs
    if (!is.null(f_fixDims)) {
      for (i in ncol(f_fixDims)) {
        f.dims[[f_fixDims[1,i]]] <- f_fixDims[2,i]
      }
    }

    # update the set of maximum dimensions
    if (!is.null(f_maxDims)) {
      for (i in ncol(f_maxDims)) {
        if (max(sp.base$f.dims[[f_maxDims[1,i]]]) == f_maxDims[2,i]) {
          f.dims[[f_maxDims[1,i]]] <- 0:f_maxDims[2,i]
        } else {
          f.dims[[f_maxDims[1,i]]] <- 1:f_maxDims[2,i]
        }
      }
    }

    # update the set of potential distances for functional inputs
    if (!is.null(f_disTypes)) {
      ids <- as.numeric(names(f_disTypes))
      for (i in length(ids)) {
        f.dist[[ids[i]]] <- f_disTypes[[i]]
      }
    }

    # update the set of potential bases for functional inputs
    if (!is.null(f_basTypes)) {
      ids <- as.numeric(names(f_basTypes))
      for (i in length(ids)) {
        f.bas[[ids[i]]] <- f_basTypes[[i]]
      }
    }

    # update the set of potential kernel functions
    if (!is.null(kerTypes)) {
      k.type <- kerTypes
    }
  }

  # fill updated space
  sp.user <- list(ds = ds, df = df, s.state = s.state, f.state = f.state,
                  f.dims = f.dims, f.dist = f.dist, f.bas = f.bas, k.type = k.type)

  return(list(sp.base = sp.base, sp.user = sp.user))
}
# ==========================================================================================================



# ==========================================================================================================
# ==========================================================================================================
getFitness <- function(model, sIn.vl = NULL, fIn.vl = NULL, sOut.vl = NULL, active = NULL) {
  # identify required statistic based on the ind.vl matrix
  if (is.null(sOut.vl)) {
    stat <- "Q2loocv"
    sOut <- model@sOut
  } else {
    stat <- "Q2hout"
    if (!is.null(active)) {
      if (length(active$s.active) > 0) sIn.pr <- sIn.vl[, active$s.active, drop = FALSE] else sIn.pr <- NULL
      if (length(active$f.active) > 0) fIn.pr <- fIn.vl[active$f.active] else fIn.pr <- NULL
    } else {
      sIn.pr <- sIn.vl
      fIn.pr <- fIn.vl
    }
    if (is.matrix(fIn.pr)) fIn.pr <- list(fIn.pr)
  }

  # compute statistic
  switch(stat,
         "Q2loocv" = {# 1: leave-one-out cross-validation Q2
           y.hat <- getOut_loocv(model)
           eta <- 1 - (mean((model@sOut - y.hat)^2)/mean((model@sOut - mean(model@sOut))^2))
         },

         "Q2hout" = {# 3: Hold-out Q2 (external validation set)
           y.hat <- quiet(predict(model, sIn.pr = sIn.pr, fIn.pr = fIn.pr)$mean)
           eta <- 1 - (mean((sOut.vl - y.hat)^2)/mean((sOut.vl - mean(sOut.vl))^2))
         })

  return(eta)
}
# ==========================================================================================================



# ==========================================================================================================
# ==========================================================================================================
splitData <- function(sIn, fIn, sOut, ind.vl) {
  ind.all <- 1:nrow(sOut) # indices of full data

  # splitting scalar inputs (if any)
  if (!is.null(sIn)) {
    sIn.tr <- sIn[ind.all[-ind.vl],,drop = FALSE]
    sIn.vl <- sIn[ind.all[ind.vl],,drop = FALSE]
  } else {
    sIn.tr <- sIn.vl <- NULL
  }

  # splitting functional inputs (if any)
  if (!is.null(fIn)) {
    fIn.tr <- lapply(fIn, function(M) M[ind.all[-ind.vl],,drop = FALSE])
    fIn.vl <- lapply(fIn, function(M) M[ind.all[ind.vl],,drop = FALSE])
  } else {
    fIn.tr <- fIn.vl <- NULL
  }

  # splitting the output
  sOut.tr <- sOut[ind.all[-ind.vl],,drop = FALSE]
  sOut.vl <- sOut[ind.all[ind.vl],,drop = FALSE]

  return(list(sIn.tr = sIn.tr, fIn.tr = fIn.tr, sOut.tr = sOut.tr,
              sIn.vl = sIn.vl, fIn.vl = fIn.vl, sOut.vl = sOut.vl))
}
# ==========================================================================================================



# ==========================================================================================================
# ==========================================================================================================
getSpacesize <- function(space) {
  # recover components
  s.state <- space$s.state
  f.state <- space$f.state
  f.dist <- space$f.dist
  f.dims <- space$f.dims
  f.bas <- space$f.bas
  k.type <- space$k.type

  # count 2 levels for each free scalar input
  n.s <- max(2 * sum(s.state == 0),1)

  # count for functional inputs
  n.fs <- rep(1, space$df)
  if (space$df > 0) {
    for (i in 1:space$df) {
      # count the number of distance types
      n.fs[i] <- n.fs[i] * length(f.dist[[i]])

      # count the number of dimensions
      n.fs[i] <- n.fs[i] * length(f.dims[[i]])

      # count the number of bases
      n.fs[i] <- n.fs[i] * length(f.bas[[i]])

      # if the input is free, add an extra level for the function inactive
      if (f.state[i] == 0) n.fs[i] <- n.fs[i] + 1
    }
  }
  n.f <- prod(n.fs)

  # count the number of kernel functions
  n.k <- length(k.type)

  # totalize and subtract 1 if all scalar and functional inputs are free
  n.solspace <- n.s * n.f * n.k - 1 * (sum(s.state, f.state) == 0)

  return(n.solspace)
}
# ==========================================================================================================



# ==========================================================================================================
# ==========================================================================================================
printSpace <- function(ds, df, space) {
  # recover components
  s.state <- space$s.state
  f.state <- space$f.state
  f.dist <- space$f.dist
  f.dims <- space$f.dims
  f.bas <- space$f.bas
  k.type <- space$k.type

  # initialize vectors of free factors and values
  free.fact <- free.values <- c()

  # initialize vectors of fixed factors and values
  fixe.fact <- fixe.values <- c()

  # check state of scalar inputs
  for (i in 1:ds) {
    if (s.state[i] == 0) {
      free.fact <- c(free.fact, paste("State of X", i, sep = ""))
      free.values <- c(free.values, "Inactive, Active")
    } else {
      fixe.fact <- c(fixe.fact, paste("State of X", i, sep = ""))
      fixe.values <- c(fixe.values, "Active")
    }
  }

  # check state of functional inputs
  for (i in 1:df) {
    if (f.state[i] == 0) {
      free.fact <- c(free.fact, paste("State of F", i, sep = ""))
      free.values <- c(free.values, "Inactive, Active")
    } else {
      fixe.fact <- c(fixe.fact, paste("State of F", i, sep = ""))
      fixe.values <- c(fixe.values, "Active")
    }
  }

  # check distances for functions
  for (i in 1:df) {
    if (length(f.dist[[i]]) > 1) {
      free.fact <- c(free.fact, paste("Dist. for F", i, sep = ""))
      free.values <- c(free.values, paste(f.dist[[i]], collapse = ", "))
    } else {
      fixe.fact <- c(fixe.fact, paste("Dist. for F", i, sep = ""))
      fixe.values <- c(fixe.values, f.dist[[i]])
    }
  }

  # check potential dimensions of functional inputs
  for (i in 1:df) {
    if (length(f.dims[[i]]) > 1) {
      free.fact <- c(free.fact, paste("Dim. of F", i, sep = ""))
      if (0 %in% f.dims[[i]]) {
        free.values <- c(free.values, paste("Orig,", paste(range(f.dims[[i]][-1]), collapse = ":")))
      } else {
        free.values <- c(free.values, paste(range(f.dims[[i]]), collapse = ":"))
      }
    } else {
      fixe.fact <- c(fixe.fact, paste("Dim. of F", i, sep = ""))
      fixe.values <- c(fixe.values, as.character(f.dims[[i]]))
    }
  }

  # check basis family for functions
  for (i in 1:df) {
    if (length(f.bas[[i]]) > 1) {
      free.fact <- c(free.fact, paste("Basis. for F", i, sep = ""))
      free.values <- c(free.values, paste(f.bas[[i]], collapse = ", "))
    } else {
      fixe.fact <- c(fixe.fact, paste("Basis. for F", i, sep = ""))
      fixe.values <- c(fixe.values, f.bas[[i]])
    }
  }

  # check kernel type
  if (length(k.type) > 1) {
    free.fact <- c(free.fact, "Kernel type")
    free.values <- c(free.values, paste(k.type, collapse = ", "))
  } else {
    fixe.fact <- c(fixe.fact, "Kernel type")
    fixe.values <- c(fixe.values, k.type)
  }

  # print table of free factors
  cat("Free")
  print(knitr::kable(cbind(free.fact, free.values), align = "c", col.names = c("Factor", "Levels"), caption = "Free factors"))

  if (length(fixe.fact) > 0) {
    # print table of fixed factors
    cat("Fixed")
    print(knitr::kable(cbind(fixe.fact, fixe.values), align = "c", col.names = c("Factor", "Levels"), caption = "Fixed factors"))
  }
  cat("\n")
}

## =============================================================================
## summary method
## =============================================================================
##' Display a summary of the structure for up to \code{n} \code{fgpm}
##' objects visited during the ACO optimization.
##'
##' @title Summary Method
##' @param object An \code{Xfgpm} object.
##' @param n Maximal number of lines (\code{fgpm} objects) to show.
##' @param ... Not used yet.
##' @return An object inheriting from \code{data.frame}.
##'
##' @method summary Xfgpm
##' @rdname summary-methods
##' @aliases summary,Xfgpm-method
##' @export
setMethod("summary", "Xfgpm",
          function(object, n = 24) {
              summary.Xfgpm(object = object, n = n)
          })

summary.Xfgpm <- function(object, n = 24, ...) {

    ds <- ncol(object@sIn)
    df <- length(object@fIn)
    n <- pmin(n, nrow(object@log.success@sols))

    ## =========================================================================
    ## When the number of variables is small enough, the structural
    ## parameters can be displayed in a single data frame. With more variables
    ## we split the content in two data frames (daf): one daf for
    ## =========================================================================

    if (4 * ds + 22 * df < 80) {
        daf <- formatShort(object@log.success@sols)
        daf <- cbind(daf,
                     "Q2" = sprintf("%5.3f", object@log.success@fitness))
        n <- pmin(n, nrow(daf))
        daf <- daf[1L:n, , drop = FALSE]
        daf <- list("Inputs and details" = daf)
    } else {
        fullDf <- object@log.success@sols
        nms <- colnames(fullDf)
        ## extract columns corresponding to the states
        indStateScal <-  grep("State_X[0-9]*", names(fullDf))
        indStateFun <-  grep("State_F[0-9]*", names(fullDf))
        indKernel <-  grep("Kernel", names(fullDf))
        activeDf <- formatShort(fullDf[ , c(indStateScal, indStateFun, indKernel),
                                       drop = FALSE])
        activeDf <- cbind(activeDf,
                          "Q2" = sprintf("%5.3f", object@log.success@fitness))
        activeDf <- activeDf[1L:n, , drop = FALSE]

        indDf <-  grep("_F[0-9]*", names(fullDf))
        funDf <- formatShort(fullDf[ , c(indDf, indKernel), drop = FALSE])
        funDf <- cbind(funDf,
                       "Q2" = sprintf("%5.3f", object@log.success@fitness))
        funDf <- funDf[1L:n, , drop = FALSE]
        daf <- list("State of inputs" = activeDf,
                 "Details for functional inputs " = funDf)
    }

    class(daf) <- c("summary.Xfgpm", "list")
    daf
}

## =============================================================================
## print method
## =============================================================================
##' @method print summary.Xfgpm
##' @export
print.summary.Xfgpm <- function(x, ...) {
    for (i in seq_along(x)) {
        cat(names(x)[i], "\n")
        print(x[[i]])
    }
}

## =============================================================================
## Re fit a fgpm model
## =============================================================================
##' Refit a \code{fgpm} model as described in a \code{Xfgpm} object.
##'
##' @title Refit a \code{fgpm} model in a \code{Xfgpm} object
##'
##' @param x A \code{Xfgpm} object.
##'
##' @param i An integer giving the index of the model to refit. The
##'     models are in decreasing fit quality as assessed by the
##'     Leave-One-Out \eqn{Q^2}{Q2}.
##'
##' @section Caution: While the syntax may suggest that the function
##'     \emph{extracts} a fitted \code{fgpm} model, this not true. The
##'     \code{fgpm} model is refitted using the call that was used
##'     when this model was assessed. The refitted \code{fgpm} model
##'     keeps the same structural parameters as the one assessed
##'     (active variables, kernel, ...), but since the optimization
##'     uses random initial values, the optimized hyper-parameters may
##'     differ from those of the corresponding \code{fgpm} in the
##'     \code{Xfgpm} object \code{x}. As a result, the model can be
##'     different and show a different LOO performance.
##'
##' @note The slot \code{@model} returns the best \code{fgpm} as
##'     assessed in a \code{Xfgm} model \code{x}. So this model can be
##'     expected to be close to the same as \code{x[[1]]}. Yet due to
##'     the refit, the two models \code{x@model} and \code{x[[1]]} can
##'     differ, see the explanations in the \bold{Caution} section.
##'
##' @export
##' @method [[ Xfgpm
##'
##' @seealso The \code{\link{modelDef}} function to extract the
##'     definition of a `fgpm` model e.g., to evaluate it using new
##'     data `sIn` `fIn`, `sOut`.
##'
##'
setMethod("[[", "Xfgpm",
          function(x, i) {
              if (i > length(x@log.success@args)) {
                  stop("'i' must be less or equal the number of 'fgpm' ",
                       "models fitted in 'x'")
              }
              ## copy the the inputs 'sIn' and 'fIn' stored in the
              ## object 'x' into the current environnement. This is
              ## necessary because these objects may have changed in
              ## the global environment.
              sIn <- x@sIn
              fIn <- x@fIn
              sOut <- x@sOut
              text <- x@log.success@args[[i]]@string
              ## text <- gsub("= sIn", "= x@sIn", text)
              ## text <- gsub("= fIn", "= x@fIn", text)
              ## cat("XXX", text, "\n")
              eval(parse(text = text)[[1]])
          })

## =============================================================================
## Retrieve a fgpm model from a Xfgpm object
## =============================================================================

##' Retrieve the \code{fgpm} model with index (or rank) \code{i} from
##' within a \code{Xfgpm} object. By evaluating this code in an
##' environment containing suitable objects \code{sIn}, \code{fIn} and
##' \code{sOut} we can re-create a \code{fgpm} object.
##'
##' The models are sorted by decreasing quality so \code{i = 1} extracts
##' the definition of the best model.
##'
##' @title Retrieve a \code{fgpm} from within a \code{Xfgpm} object
##'
##' @param object A \code{Xfgpm} object as created by
##' \code{\link{fgpm_factory}}.
##'
##' @param ind The index (or rank) of the model in \code{object}.
##'
##' @return A parsed R code defining the \code{fgpm} model.
##'
##' @note Remind that the models are sorted by decreasing quality so
##'     \code{i = 1} extracts the definition of the best model.
##'
##' @seealso the \code{\link{[[,Xfgpm-method}} that can also be used
##'     to re-create a \code{fgpm} object using \emph{the same data}
##'     as that used to create the \code{Xfgpm} object in
##'     \code{object}.
##'
##' @export
##'
##' @examples
##' set.seed(100)
##' n.tr <- 32
##' x1 <- x2 <- x3 <- x4 <- x5 <- seq(0, 1, length = n.tr^(1/5))
##'
##' sIn <- as.matrix(expand.grid(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5))
##' fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
##'             f2 = matrix(runif(n.tr * 22), ncol = 22))
##' sOut <- fgp_BB7(sIn, fIn, n.tr)
##' xm <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut)
##'
##' ## 'xm@model' is the best 'fgpm' model in 'xm'
##' plot(xm@model)
##' modelDef(xm, i = 1)
##'
##' ## Define new data
##' n.new <- 3^5
##' x1 <- x2 <- x3 <- x4 <- x5 <- seq(0, 1, length = n.new^(1/5))
##'
##' ## replace the data objects from which the model is created
##' sIn <- as.matrix(expand.grid(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5))
##' fIn <- list(f1 = matrix(runif(n.new * 10), ncol = 10),
##'             f2 = matrix(runif(n.new * 22), ncol = 22))
##' sOut <- fgp_BB7(sIn, fIn, n.new)
##' ## Now evaluate the definition. Since the objects 'sIn' and 'fIn'
##' ## have been replaced by new values, the perfomance will differ.
##' fgpm.new <- eval(modelDef(xm, i = 1))
##' plot(fgpm.new, main = "Re-created 'fgpm' model with different data")
##' plot(xm[[1]], main = "Re-created 'fgpm' model with the same data")
modelDef <- function(object, ind) {
    parse(text = object@log.success@args[[ind]]@string[[1]])
}
