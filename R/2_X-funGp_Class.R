#' @title Class: data structures related to the kernel of a funGp model
#' @description Fill this!!!!!!!!!
#'
#' @slot model Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#' @slot fitness Hola.
#' @slot structure Hola.
#' @slot factoryCall Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#' @slot method Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#' @slot stat Hola.
#' @slot details Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#' @slot log Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#'
#' @include 2_funGp_Class.R
#' @include 3_ant_search.R
#' @include 8_outilsCode.R
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setClass("X-funGp",
         representation(
           model = "funGp",                # kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}
           fitness = "numeric",            # model fitness
           structure = "data.frame",       # model fitness
           factoryCall = "factoryCall",    # distance type. To be chosen from {"scalar", "functional"}
           method = "character",           # search method
           stat = "character",             # search method
           details = "list",               # search method
           log = "antsLog"                 # search method
         ),
         validity = function(object) {T})


# funGp master function: used for construction and training of a funGP model
# ----------------------------------------------------------------------------------------------------------
#' @title Fill!!!
#' @description Fill. This.
#' @param sIn fill
#' @param fIn fill
#' @param sOut fill
#' @param constraints fill
#' @param setup fill
#' @param nugget Fill!!!!!!!!!!
#' @param n.starts Fill!!!!!!!!!!
#' @param n.presample Fill!!!!!!!!!!
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
#' @keywords internal
funGp_factory <- function(sIn = NULL, fIn = NULL, sOut = NULL, ind.vl = NULL,
                          method = "ACO", constraints = list(), setup = list(),
                          nugget = 1e-8, n.starts = 1, n.presample = 20) {
  # user inputs, remove this at the end !!!!!!!!!!!!!!!!!!!!!!!!!
  s_keepActive <- c(1) # keep X1 always active
  f_keepActive <- NULL # Do not keep any functional input always active, let them free
  f_fixDims <- matrix(c(2,4), ncol = 1) # set dimension of F2 at 4
  f_maxDims <- matrix(c(1,5), ncol = 1) # set max dimension of F1 at 5
  f_disTypes <- list("2" = c("L2_byindex")) # only test L2_index distance for F2
  f_basTypes <- list("1" = c("B-splines")) # only B-splines projection for F1
  kerTypes <- c("matern5_2", "matern3_2") # test only matern kernels
  constraints <- list(s_keepActive = s_keepActive, f_keepActive = f_keepActive, f_fixDims = f_fixDims, f_maxDims = f_maxDims,
                      f_disTypes = f_disTypes, f_basTypes = f_basTypes, kerTypes = kerTypes)

  # extra arguments for model call
  extargs <- list(nugget = nugget, n.starts = n.starts, n.presample = n.presample)

  # define solution space based on user inputs
  solspace <- setSpace(sIn, fIn, constraints)

  if (all(!is.null(sIn), !is.null(fIn))) { # Hybrid-input case *******************************************
    print("I'm hybrid!")
  } else if(!is.null(fIn)) { # functional-input case ***************************************
    print("I'm functional!")
  } else if(!is.null(sIn)) { # scalar-input case *******************************************
    print("I'm scalar!")
  } else { # error: no inputs were provided
    stop("The user must provide either a scalar-input matrix, a functional-input list or both of them. None has been detected.")
  }

  # prepare input and output structures based on case: 1. LOOCV; 2. HOUT
  # if (is.null(ind.vl)) { # case 1
  #   sIn.tr <- sIn
  #   fIn.tr <- fIn
  #   sOut.tr <- sOut
  # } else { # case 2
  # ind.vl <- as.matrix(ind.vl)
  #   sIn.tr <- sIn
  #   fIn.tr <- fIn
  #   sOut.tr <- sOut
  # }

  if (!is.null(ind.vl)) {
    ind.vl <- as.matrix(ind.vl)
    stat <- paste("Q2hout.", (nrow(sOut) - nrow(ind.vl)), ".", nrow(ind.vl), ".", ncol(ind.vl), sep = "")
  } else {
    stat <- "Q2loocv"
  }


  # optimize model structure
  switch(method,
         "ACO" = {# 1: Ant Colony Optimization
           opt <- master_ACO(sIn, fIn, sOut, ind.vl, solspace, setup, extargs)
         },

         "ES" = {# 3: Exhaustive Search
           # opt <- matern32_cor(Ms, thetas)
         })

  X.model <- new("X-funGp")
  X.model@model <- opt$model
  X.model@fitness <- opt$b.fitness
  X.model@structure <- opt$sol.vec
  X.model@factoryCall@string <- gsub("^ *|(?<= ) | *$", "", paste0(deparse(match.call()), collapse = " "), perl = T)
  X.model@method <- "Ants"
  X.model@stat <- stat
  X.model@log <- opt$log
  X.model@details <- opt$details

  return(X.model)
}

setSpace <- function(sIn, fIn, constraints) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  # recover individual constraints from user inputs
  s_keepActive <- constraints$s_keepActive
  f_keepActive <- constraints$f_keepActive
  f_fixDims <- constraints$f_fixDims
  f_maxDims <- constraints$f_maxDims
  f_disTypes <- constraints$f_disTypes
  f_basTypes <- constraints$f_basTypes
  kerTypes <- constraints$kerTypes

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
  k.type <- c("gaussian", "matern5_2", "matern3_2")
  sp.base <- list(s.state = s.state, f.state = f.state, f.dims = f.dims, f.dist = f.dist, f.bas = f.bas, k.type = k.type)

  # modify the solution space if the user has specified any constraint
  if (length(constraints) > 0) {
    # update state of scalar inputs
    if (!is.null(s_keepActive)) s.state[s_keepActive] <- 1

    # update state of functional inputs
    if (!is.null(f_keepActive)) f.state[f_keepActive] <- 1

    # update the set of potential dimensions for functional inputs
    if (!is.null(f_fixDims)) {
      for (i in ncol(f_fixDims)) {
        f.dims[[f_fixDims[1,i]]] <- f_fixDims[2,i]
      }
    }

    # update the set of maximum dimensions
    if (!is.null(f_maxDims)) {
      for (i in ncol(f_maxDims)) {
        f.dims[[f_maxDims[1,i]]] <- 1:f_maxDims[2,i]
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

getFitness <- function(model, sIn.vl = NULL, fIn.vl = NULL, sOut.vl = NULL, active = NULL) {
  # identify required statistic based on the ind.vl matrix
  if (is.null(sOut.vl)) {
    stat <- "Q2loocv"
    sOut <- model@sOut
  } else {
    stat <- "Q2hout"
    if (length(active$s.active) > 0) sIn.pr <- sIn.vl[, active$s.active, drop = F] else sIn.pr <- NULL
    if (length(active$f.active) > 0) fIn.pr <- fIn.vl[active$f.active] else fIn.pr <- NULL
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

splitData <- function(sIn, fIn, sOut, ind.vl) {
  ind.all <- 1:nrow(sOut) # indices of full data

  # splitting scalar inputs (if any)
  if (!is.null(sIn)) {
    sIn.tr <- sIn[ind.all[-ind.vl],,drop = F]
    sIn.vl <- sIn[ind.all[ind.vl],,drop = F]
  } else {
    sIn.tr <- sIn.vl <- NULL
  }

  # splitting functional inputs (if any)
  if (!is.null(fIn)) {
    fIn.tr <- lapply(fIn, function(M) M[ind.all[-ind.vl],,drop = F])
    fIn.vl <- lapply(fIn, function(M) M[ind.all[ind.vl],,drop = F])
  } else {
    fIn.tr <- fIn.vl <- NULL
  }

  # splitting the output
  sOut.tr <- sOut[ind.all[-ind.vl],,drop = F]
  sOut.vl <- sOut[ind.all[ind.vl],,drop = F]

  return(list(sIn.tr = sIn.tr, fIn.tr = fIn.tr, sOut.tr = sOut.tr,
              sIn.vl = sIn.vl, fIn.vl = fIn.vl, sOut.vl = sOut.vl))
}

printSpace <- function(ds, df, space) {
  # recover components
  s.state <- space$s.state
  f.state <- space$f.state
  f.dims <- space$f.dims
  f.dist <- space$f.dist
  f.fam <- space$f.fam
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

  # check basis family for functions
  for (i in 1:df) {
    if (length(f.fam[[i]]) > 1) {
      free.fact <- c(free.fact, paste("Basis. for F", i, sep = ""))
      free.values <- c(free.values, paste(f.fam[[i]], collapse = ", "))
    } else {
      fixe.fact <- c(fixe.fact, paste("Basis. for F", i, sep = ""))
      fixe.values <- c(fixe.values, f.fam[[i]])
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
