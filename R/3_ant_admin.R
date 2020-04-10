# ==========================================================================================================
# S4 class for the log of fgpm_factory optimizations using the ACO algorithm
# ==========================================================================================================
#' @title S4 class for log of models explored by ant colony in funGp
#' @description Register of model structures and their performance statistic, if available.
#'
#' @slot sols Object of class \code{"data.frame"}. Compendium of model structures arranged by rows. Each
#'   column is linked to one structural parameter of the model such as the state of one variable (inactive,
#'   active) or the type of kernel function.
#' @slot args Object of class \code{"list"}. Compendium of model structures represented by objects of class
#'   \code{"\linkS4class{modelCall}"}
#' @slot fitness Object of class \code{"numeric"}. Performance statistic of each model, if available.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @rdname antsLog-class
#' @export
setClass("antsLog",
         representation(
           sols = "data.frame",            # compendium of model structures with selections as factor levels
           args = "list",                  # compendium of model structures with selections as modelCall objects
           fitness = "numeric"             # performance statistic of each model, if available
         ),
         validity = function(object) {TRUE})
# ==========================================================================================================



# ==========================================================================================================
# Master function to process ACO tasks
# ==========================================================================================================
master_ACO <- function(sIn, fIn, sOut, ind.vl, solspace, setup, extargs, time.str, time.lim, quietly, par.clust) {
  # set heuristic parameters based on defaults and user specifications
  param <- setParams_ACO(setup, length(fIn))

  # set pheromones based on the solution space and initial parameters
  phero <- setEnvir_ACO(solspace, param)

  # perform exploration
  res <- run_ACO(sIn, fIn, sOut, ind.vl, param, phero, solspace$sp.base, extargs, time.str, time.lim, quietly, par.clust)

  # correct the howCalled of best model
  res$model@howCalled@string <- getCall_ACO(res$sol.vec, sIn, fIn, res$sol.args, solspace$sp.base, extargs)

  # format solution vector
  res$sol.vec <- vec2DFrame_ACO(res$sol.vec, sIn, fIn)

  # fill antLog for successfull ants
  res$log.suc <- getLog_ACO(sIn, fIn, res$log.suc, res$log.fitness, solspace$sp.base, extargs)

  # fill antLog for crashed ants if any
  if(length(res$log.cra) > 0) {
    res$log.cra <- getLog_ACO(sIn, fIn, res$log.cra, numeric(), solspace$sp.base, extargs)
  } else {
    res$log.cra <- new("antsLog")
  }

  return(res)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to prepare the parameters of the heuristic considering the user inputs
# ==========================================================================================================
setParams_ACO <- function(setup, df) {
  if (!is.null(setup$n.iter)) n.iter <- setup$n.iter else n.iter <- 15
  if (!is.null(setup$n.pop)) n.pop <- setup$n.pop else n.pop <- 10
  if (!is.null(setup$tao0)) tao0 <- setup$tao0 else tao0 <- .1
  if (!is.null(setup$dop.s)) dop.s <- setup$dop.s else dop.s <- 1#tao0
  if (!is.null(setup$dop.f)) dop.f <- setup$dop.f else dop.f <- 1#tao0
  if (!is.null(setup$delta.f)) {
    if (length(setup$delta.f) == df) {
      delta.f <- setup$delta.f
    } else if (length(setup$delta.f) == 1) {
      delta.f <- rep(setup$delta.f, df)
    } else {
      stop("The parameter delta.f should be either of length one or of length equal to the number of functional inputs.")
    }
  } else {
    delta.f <- rep(2, max(1,df))
  }
  if (!is.null(setup$dispr.f)) {
    if (length(setup$dispr.f) == df) {
      dispr.f <- setup$dispr.f
    } else if (length(setup$dispr.f) == 1) {
      dispr.f <- rep(setup$dispr.f, df)
    } else {
      stop("The parameter dispr.f should be either of length one or of length equal to the number of functional inputs.")
    }
  } else {
    dispr.f <- rep(1.4, max(1,df))
  }
  if (!is.null(setup$q0)) q0 <- setup$q0 else q0 <- .95
  if (!is.null(setup$rho.l)) rho.l <- setup$rho.l else rho.l <- .1
  if (!is.null(setup$u.gbest)) u.gbest <- setup$u.gbest else u.gbest <- F
  if (!is.null(setup$n.ibest)) n.ibest <- setup$n.ibest else n.ibest <- 1
  if (!is.null(setup$rho.g)) rho.g <- setup$rho.g else rho.g <- .1

  params <- list(n.iter = n.iter, n.pop = n.pop, tao0 = tao0, dop.s = dop.s, dop.f = dop.f,
                 delta.f = delta.f, dispr.f = dispr.f, q0 = q0, rho.l = rho.l,
                 u.gbest = u.gbest, n.ibest = n.ibest, rho.g = rho.g)

  return(params)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to generate the initial pheromone load
# ==========================================================================================================
#' @importFrom scales rescale
setEnvir_ACO <- function(solspace, param) {
  # recover base and user-defined solution space
  sp.base <- solspace$sp.base
  sp.user <-solspace$sp.user

  # recover base solution space
  s.state.0 <- sp.base$s.state
  f.state.0 <- sp.base$f.state
  f.dims.0 <- sp.base$f.dims
  f.dist.0 <- sp.base$f.dist
  f.bas.0 <- sp.base$f.bas
  k.type.0 <- sp.base$k.type

  # recover user-defined solution space
  s.state.u <- sp.user$s.state
  f.state.u <- sp.user$f.state
  f.dims.u <- sp.user$f.dims
  f.dist.u <- sp.user$f.dist
  f.bas.u <- sp.user$f.bas
  k.type.u <- sp.user$k.type

  # recover required parameters from heuristic setup
  tao0 <- param$tao0
  dop.s <- param$dop.s
  dop.f <- param$dop.f
  delta.f <- param$delta.f
  dispr.f <- param$dispr.f

  # extract input dimensions
  ds <- sp.user$ds
  df <- sp.user$df

  # create pheromone list
  phero <- list()
  layers <- c()
  if (ds > 0) {
    # set up pheromones related to scalar inputs
    if (s.state.u[1] == 0) {
      phero[[1]] <- matrix(c(tao0, tao0/dop.s), ncol = 2, byrow = TRUE)
    } else {
      phero[[1]] <- matrix(c(1, 0), ncol = 2) # set up initial pheromones
    }
    rownames(phero[[1]]) <- "Orig"
    colnames(phero[[1]]) <- c("Active", "Inactive")
    layers <- c(layers, paste("State X", 1, sep =""))

    if (ds > 1) {
      for (i in 2:ds){
        # state of scalar input i
        # ____________________________________________________________________________________
        if (s.state.u[i] == 0) {
          phero[[i]] <- matrix(rep(c(tao0, tao0/dop.s), 2), nrow = 2, byrow = TRUE)
        } else {
          phero[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = TRUE) # set up initial pheromones
        }
        # set up colnames
        colnames(phero[[i]]) <- rownames(phero[[i]]) <- c("Active", "Inactive")
        # ____________________________________________________________________________________
        layers <- c(layers, paste("State X", i, sep = ""))
      }
    }
  }

  # set up visibility and pheromones related to functional inputs
  if (df > 0) {
    for (j in 1:df) {
      # state of functional input j
      # ____________________________________________________________________________________
      i <- length(phero) + 1
      if (all(i == 1, j == 1)) {
        if (f.state.u[j] == 0) { # free, distribute pheromones and assign default visibilities
          phero[[i]] <- matrix(c(tao0, tao0/dop.f), ncol = 2)
        } else { # fixed, put all load in first colum
          phero[[i]] <- matrix(c(1, 0), ncol = 2) # set up initial pheromones
        }
        rownames(phero[[i]]) <- "Orig"
        colnames(phero[[i]]) <- c("Active", "Inactive")

      } else {
        if (f.state.u[j] == 0) { # free, distribute pheromones and assign default visibilities
          phero[[i]] <- matrix(rep(c(tao0, tao0/dop.f), 2), nrow = 2, byrow = TRUE)
        } else { # fixed, put all load in first colum
          phero[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = TRUE)
        }
        # <--> names
        rownames(phero[[i]]) <- colnames(phero[[i-1]])
        colnames(phero[[i]]) <- c("Active", "Inactive")
      }
      layers <- c(layers, paste("State F", j, sep =""))
      # ____________________________________________________________________________________


      # distance for functional input j
      # ____________________________________________________________________________________
      i <- length(phero) + 1 # i <- length(visib) + 1
      nr <- ncol(phero[[i-1]]) # nr <- ncol(visib[[i-1]])
      # <--> identify active levels
      v <- rep(0, length(f.dist.0[[j]]))
      t <- unname(sapply(f.dist.u[[j]], function(x) which(f.dist.0[[j]] == x)))
      # <--> assign initial pheromones
      if (length(t) == 1)  v[t] <- 1 else v[t] <- tao0
      phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = TRUE)
      # <--> names
      rownames(phero[[i]]) <- colnames(phero[[i-1]])
      colnames(phero[[i]]) <- f.dist.0[[j]]
      layers <- c(layers, paste("Dist. FI", j, sep =""))
      # ____________________________________________________________________________________


      # dimension for functional input j
      # ____________________________________________________________________________________
      i <- length(phero) + 1 # i <- length(visib) + 1
      nr <- ncol(phero[[i-1]]) # nr <- ncol(visib[[i-1]])
      # <--> identify active levels
      vgro <- vind <- rep(0, length(f.dims.0[[j]])) # zero vector with as many elements as potential dimensions in the original set
      t <- unname(sapply(f.dims.u[[j]], function(x) which(f.dims.0[[j]] == x))) # locations of dimensions still active
      # <--> assign initial pheromones
      if (length(t) == 1) {
        vind[t] <- vgro[t] <- 1
      } else {
        vgro[t] <- tao0
        if (0 %in% f.dims.u[[j]]) {
          l <- decay(tao0 = tao0, delta = delta.f[j], dispr = dispr.f[j], k = max(f.dims.u[[j]]), doplot = FALSE, deliver = TRUE)
        } else {
          l <- decay(tao0 = tao0, delta = delta.f[j], dispr = dispr.f[j], k = max(f.dims.u[[j]]), doplot = FALSE, deliver = TRUE)[-1]
        }
        vind[t] <- l
      }
      # <--> wrap visibility
      phero[[i]] <- rbind(vgro, vind)
      # <--> names
      rownames(phero[[i]]) <- colnames(phero[[i-1]])
      colnames(phero[[i]]) <- f.dims.0[[j]]
      layers <- c(layers, paste("Dim. FI", j, sep =""))
      # ____________________________________________________________________________________


      # basis for functional input j
      # ____________________________________________________________________________________
      i <- length(phero) + 1
      nr <- ncol(phero[[i-1]])
      # <--> identify active levels
      v <- rep(0, length(f.bas.0[[j]])) # zero vector with as many elements as basis families in the original set
      t <- unname(sapply(f.bas.u[[j]], function(x) which(f.bas.0[[j]] == x))) # locations of basis still active
      # <--> assign initial pheromones
      if (length(t) == 1)  v[t] <- 1 else v[t] <- tao0
      phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = TRUE)
      # <--> names
      colnames(phero[[i]]) <- f.bas.0[[j]]
      rownames(phero[[i]]) <- f.dims.0[[j]]
      layers <- c(layers, paste("Basis FI", j, sep =""))
      # ____________________________________________________________________________________
    }
  }

  # kernel function
  # ____________________________________________________________________________________
  i <- length(phero) + 1
  nr <- ncol(phero[[i-1]])
  # <--> identify active levels
  v <- rep(0, length(k.type.0)) # zero vector with as many elements as kernel functions in the original set
  t <- unname(sapply(k.type.u, function(x) which(k.type.0 == x))) # locations of kernels still active
  # <--> assign initial pheromones
  if (length(t) == 1)  v[t] <- 1 else v[t] <- tao0
  phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = TRUE)
  # <--> names
  rownames(phero[[i]]) <- colnames(phero[[i-1]])
  colnames(phero[[i]]) <- k.type.0
  layers <- c(layers, "Kernel")
  # ____________________________________________________________________________________

  # set names to pheromone list
  names(phero) <- layers

  # return(list(phero = phero, visib = visib))
  return(phero)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to transform ants into data structures for fgpm
# ==========================================================================================================
formatSol_ACO <- function(ant, sIn, fIn, base) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  # recover base components
  s.state <- base$s.state
  f.state <- base$f.state
  f.dims <- base$f.dims
  f.dist <- base$f.dist
  f.bas <- base$f.bas
  k.type <- base$k.type

  # remove inactive scalar variables
  if (ds > 0) {
    s.active.ac <- which(ant[1:ds] == 1) # index of active scalar inputs
    if (length(s.active.ac) > 0) {
      sIn.ac <- sIn[,s.active.ac, drop = FALSE]
    } else {
      sIn.ac <- NULL
    }
  } else {
    sIn.ac <- NULL
  }
  piv <- ds

  defargs <- formals(fgpm) # in case they are required below
  if (df > 0) {
    f.active.ac <- c() # index of active functional inputs
    f.dist.ac <- rep(0, df) # distance for each active functional input
    f.dims.ac <- rep(0, df) # dimension of each active functional input
    f.bas.ac <- rep(0, df) # basis function for each active functional input
    for (j in 1:df) {
      piv <- piv + 1
      if (ant[piv] == 1) f.active.ac <- c(f.active.ac, j)
      piv <- piv + 1
      if (ant[piv] == -1) { # assign default value
        f.dist.ac[[j]] <- defargs$f_disType
        f.dims.ac[[j]] <- defargs$f_pdims
        f.bas.ac[[j]] <- defargs$f_basType
        piv <- piv + 2
      } else {
        f.dist.ac[[j]] <- f.dist[[j]][ant[piv]]
        piv <- piv + 1
        f.dims.ac[[j]] <- ant[piv]
        piv <- piv + 1
        f.bas.ac[[j]] <- f.bas[[j]][ant[piv]]
      }
    }
    if (length(f.active.ac) > 0) {
      fIn.ac <- fIn[f.active.ac]
      f.dist.ac <- f.dist.ac[f.active.ac]
      f.dims.ac <- f.dims.ac[f.active.ac]
      f.bas.ac <- f.bas.ac[f.active.ac]

    } else { # set by default values
      fIn.ac <- defargs$fIn
      f.dist.ac <- defargs$f_disType
      f.dims.ac <- defargs$f_pdims
      f.bas.ac <- defargs$f_basType
    }
  } else { # set by default values
    fIn.ac <- defargs$fIn
    f.dist.ac <- defargs$f_disType
    f.dims.ac <- defargs$f_pdims
    f.bas.ac <- defargs$f_basType
    piv <- piv + df * 4
  }
  piv <- piv + 1
  k.type.ac <- k.type[ant[piv]]

  return(list(sIn = sIn.ac, fIn = fIn.ac, kerType = k.type.ac, f_disType = f.dist.ac, f_pdims = f.dims.ac, f_basType = f.bas.ac))
}
# ==========================================================================================================



# ==========================================================================================================
# Function to identify active functional inputs in an ant
# ==========================================================================================================
getActiveIn_ACO <- function(ant, sIn, fIn, base) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  # recover base components
  s.state <- base$s.state
  f.state <- base$f.state

  # index of active scalar inputs
  s.active.ac <- NULL
  if (ds > 0) {
    act <- unname(which(ant[1:ds] == 1))
    if (length(act) > 0) {
      s.active.ac <- act
    }
  }

  # index of active functional inputs
  f.active.ac <- NULL
  if (df > 0) {
    act <- unname(which(ant[grepl("State F", names(ant))] == 1))
    if (length(act) > 0) {
      f.active.ac <- act
    }
  }

  return(list(s.active = s.active.ac, f.active = f.active.ac))
}
# ==========================================================================================================



# ==========================================================================================================
# Function to reconstructu the fgpm call for a given ant
# ==========================================================================================================
getCall_ACO <- function(ant, sIn, fIn, args, base, extargs) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  # identify used inputs
  in.used <- getActiveIn_ACO(ant, sIn, fIn, base)
  s.used <- in.used$s.active
  f.used <- in.used$f.active

  # set string for scalar inputs
  if (!is.null(s.used)) {
    if (length(s.used) == ncol(sIn)) {
      s.str <- "sIn = sIn"
    } else if (length(s.used) == 1) {
      s.str <- paste("sIn = sIn[,", s.used, ",drop=F]", sep = "")
    } else if (all(diff(s.used) == 1)) {
      s.str <- paste("sIn = sIn[,", paste(range(s.used), collapse = ":"), "]", sep = "")
    } else {
      s.str <- paste("sIn = sIn[,c(", paste(s.used, collapse = ","), ")]", sep = "")
    }
  } else {
    s.str <- NULL
  }

  # set string for functional inputs
  if (!is.null(f.used)) {
    if (length(f.used) == length(fIn)) {
      f.str <- "fIn = fIn"
    } else if (length(f.used) == 1) {
      f.str <- paste("fIn = fIn[", f.used, "]", sep = "")
    } else if (all(diff(f.used) == 1)) {
      f.str <- paste("fIn = fIn[", paste(range(f.used), collapse = ":"), "]", sep = "")
    } else {
      f.str <- paste("fIn = fIn[c(", paste(f.used, collapse = ","), ")]", sep = "")
    }

    # prepare f_disType string
    if (length(args$f_disType) == 1) {
      f_disType.str <- paste('f_disType = "', args$f_disType, '"', sep = "")
    } else if (length(unique(args$f_disType)) > 1) {
      f_disType.str <- paste('f_disType = c("', paste(args$f_disType, collapse = '", "'), '")', sep = "")
    } else {
      f_disType.str <- paste('f_disType = "', args$f_disType[1], '"', sep = "")
    }

    # prepare f_pdims string
    if (length(args$f_pdims) == 1) {
      f_pdims.str <- paste("f_pdims = ", args$f_pdims, sep = "")
    } else if (length(unique(args$f_pdims)) > 1) {
      f_pdims.str <- paste('f_pdims = c(', paste(args$f_pdims, collapse = ', '), ')', sep = "")
    } else {
      f_pdims.str <- paste("f_pdims = ", args$f_pdims[1], sep = "")
    }

    # prepare f_basType string
    if (length(args$f_basType) == 1) {
      f_basType.str <- paste('f_basType = "', args$f_basType, '"', sep = "")
    } else if (length(unique(args$f_basType)) > 1) {
      f_basType.str <- paste('f_basType = c("', paste(args$f_basType, collapse = '", "'), '")', sep = "")
    } else {
      f_basType.str <- paste('f_basType = "', args$f_basType[1], '"', sep = "")
    }

  } else {
    f.str <- f_disType.str <- f_pdims.str <- f_basType.str <- NULL
  }

  # prepare sOut string
  sy.str <- "sOut = sOut"

  # prepare kerType string
  kerType.str <- paste('kerType = "', args$kerType, '"', sep = "")

  # set strings for extra arguments: nugget, n.starts, n.presample
  defargs <- formals(fgpm)
  if (extargs$nugget != defargs$nugget) nugg.str <- paste("nugget = ", extargs$nugget, sep = "") else nugg.str <- NULL
  if (extargs$n.starts != defargs$n.starts) n.starts.str <- paste("n.starts = ", extargs$n.starts, sep = "") else n.starts.str <- NULL
  if (extargs$n.presample != defargs$n.presample) {
    n.presample.str <- paste("n.presample = ", extargs$n.presample, sep = "")
  } else {
    n.presample.str <- NULL
  }

  # merge strings to produce model call
  full.str <- c(s.str, f.str, sy.str, f_disType.str, f_pdims.str, f_basType.str, kerType.str, nugg.str, n.starts.str, n.presample.str)
  modcall <- paste("fgpm(", paste(full.str, collapse = ", "), ")", sep = "")

  return(modcall)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to translate ants coded levels to lexicographic tags in data.frame format
# ==========================================================================================================
vec2DFrame_ACO <- function(sol.vec, sIn, fIn) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  names <- vals <- rep(0,length(sol.vec))
  if (ds > 0) {
    for (i in 1:ds) {
      names[i] <- paste("State_X", i, sep = "")
      vals[i] <- c("On", "Off")[sol.vec[i]]
    }
    piv <- ds
  }

  if (df > 0) {
    piv <- max(0, ds)
    for (i in 1:df) {
      piv <- piv + 1
      names[piv] <- paste("State_F", i, sep = "")
      vals[piv] <- c("On", "Off")[sol.vec[piv]]

      piv <- piv + 1
      names[piv] <- paste("Distance_F", i, sep = "")
      if (sol.vec[piv] < 0) vals[piv] <- "--" else vals[piv] <- c("L2_bygroup", "L2_byindex")[sol.vec[piv]]

      piv <- piv + 1
      names[piv] <- paste("Dim_F", i, sep = "")
      if (sol.vec[piv] < 0) vals[piv] <- "-" else vals[piv] <- sol.vec[piv]

      piv <- piv + 1
      names[piv] <- paste("Prj_basis_F", i, sep = "")
      if (sol.vec[piv] < 0) vals[piv] <- "--" else vals[piv] <- c("B-splines", "PCA")[sol.vec[piv]]
    }
  }

  piv <- piv + 1
  names[piv] <- "Kernel"
  vals[piv] <- c("gauss", "matern5_2", "matern3_2")[sol.vec[piv]]

  # data.frame-ize
  vals <- data.frame(matrix(vals, nrow = 1))
  names(vals) <- names

  return(vals)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to build ant logs (success and crashes)
# ==========================================================================================================
getLog_ACO <- function(sIn, fIn, log.vec, log.fitness, base, extargs) {
  mylog <- new("antsLog")
  args <- list()
  for (i in 1:nrow(log.vec)) {
    mc <- new("modelCall")
    mc@string <- getCall_ACO(log.vec[i,], sIn, fIn, formatSol_ACO(log.vec[i,], sIn, fIn, base), base, extargs)
    args[[i]] <- mc
  }
  mylog@args <- args
  mylog@sols <- do.call("rbind", apply(log.vec, 1, function(x) vec2DFrame_ACO(x, sIn, fIn)))
  mylog@fitness <- log.fitness
  return(mylog)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to evaluate ants for the Q2loocv statistic
# ==========================================================================================================
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan cluster
#' @importFrom progressr with_progress progressor
eval_loocv_ACO <- function(sIn, fIn, sOut, extargs, base, ants, time.str, time.lim, quietly, par.clust) {
  # recover population size
  n.pop <- nrow(ants)

  if (!is.null(par.clust)) {
    # register parallel backend
    registerDoFuture()
    plan(cluster, workers = par.clust)

    # evaluate the ants as models
    with_progress({
      p <- progressor(along = 1:n.pop, auto_finish = FALSE)
      result <- foreach(i = 1:n.pop, .errorhandling = "pass") %dopar% {
        dt <- difftime(Sys.time(), time.str, units = 'secs')
        if (dt < time.lim) {
          if (quietly) {
            amf <- quiet(fitNtests_ACO(ants[i,], sIn, fIn, sOut, extargs, base))
          } else {
            amf <- fitNtests_ACO(ants[i,], sIn, fIn, sOut, extargs, base)
          }
          p()
          return(amf)
        } else {
          p()
          return(NULL)
        }
      }
    })

  } else {
    # set up progress bar
    pb <- txtProgressBar(min = 0, max = n.pop, style = 3)

    # evaluate the ants as models
    result <- list()
    for (i in 1:n.pop) {
      if (quietly) {
        result[[i]] <- quiet(fitNtests_ACO(ants[i,], sIn, fIn, sOut, extargs, base))
      } else {
        cat("\n")
        result[[i]] <- fitNtests_ACO(ants[i,], sIn, fIn, sOut, extargs, base)
      }

      setTxtProgressBar(pb, i)

      # check if we are still on time
      dt <- difftime(Sys.time(), time.str, units = 'secs')
      if (dt >= time.lim) {
        if (i < n.pop) {
          result[(i+1):n.pop] <- NULL
          setTxtProgressBar(pb, n.pop)
        }
        break
      }
    }
    close(pb)
  }
  return(result)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to evaluate ants for the Q2hout statistic
# ==========================================================================================================
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan cluster
#' @importFrom progressr with_progress progressor
eval_houtv_ACO <- function(sIn, fIn, sOut, extargs, base, ants, ind.vl, time.str, time.lim, quietly, par.clust) {
  # recover population size and number of replicates
  n.pop <- nrow(ants)
  n.rep <- ncol(ind.vl)

  if (!is.null(par.clust)) {
    # register parallel backend
    registerDoFuture()
    plan(cluster, workers = par.clust)

    # evaluate the ants as models
    with_progress({
      p <- progressor(along = 1:n.pop, auto_finish = FALSE)
      result <- foreach(i = 1:n.pop, .errorhandling = "pass") %dopar% {
        dt <- difftime(Sys.time(), time.str, units = 'secs')
        if (dt < time.lim) {
          sub_result <- list()
          for (j in 1:n.rep) {
            # split data into training and validation
            data <- splitData(sIn, fIn, sOut, ind.vl[,j])

            # identify active inputs of both types
            active <- getActiveIn_ACO(ants[i,], sIn, fIn, base)

            if (quietly) {
              sub_result[[j]] <- quiet(fitNtests_ACO(ants[i,], data$sIn.tr, data$fIn.tr, data$sOut.tr, extargs, base,
                                                     ind.vl, data$sIn.vl, data$fIn.vl, data$sOut.vl, active))
            } else {
              cat("\n")
              sub_result[[j]] <- fitNtests_ACO(ants[i,], data$sIn.tr, data$fIn.tr, data$sOut.tr, extargs, base,
                                               ind.vl, data$sIn.vl, data$fIn.vl, data$sOut.vl, active)
            }

            # check if we are still on time
            dt <- difftime(Sys.time(), time.str, units = 'secs')
            if (dt >= time.lim) {
              if (j < n.rep) {
                sub_result[(j+1):n.rep] <- NULL
              }
              break
            }
          }

          # extract complete evaluations
          done <- which(!sapply(sub_result, is.null))
          argsList <- lapply(sub_result[done], `[[`, 1)
          modelList <- lapply(sub_result[done], `[[`, 2)
          fitnessVec <- sapply(sub_result[done], `[[`, 3)

          # identify crashes and usable models
          ids.cr <- which(is.na(fitnessVec))
          ids.ok <- which(!is.na(fitnessVec))

          # remove crashes and compute average model fitness
          if (length(ids.ok) > 0) {
            fitness <- mean(fitnessVec[ids.ok])
            piv <- sample.vec(which(fitnessVec == max(fitnessVec)))
            args <- argsList[[piv]]
            model <- modelList[[piv]]
            amf <- list(args, model, fitness)
          } else {
            amf <- NULL
          }
          p()
          return(amf)
        } else {
          p()
          return(NULL)
        }
      }
    })

  } else {
    # set up progress bar
    pb <- txtProgressBar(min = 0, max = n.pop, style = 3)

    # evaluate the ants as models
    result <- list()
    for (i in 1:n.pop) {
      sub_result <- list()
      for (j in 1:n.rep) {
        # split data into training and validation
        data <- splitData(sIn, fIn, sOut, ind.vl[,j])
        # identify active inputs of both types
        active <- getActiveIn_ACO(ants[i,], sIn, fIn, base)
        # evaluate model
        if (quietly) {
          sub_result[[j]] <- quiet(fitNtests_ACO(ants[i,], data$sIn.tr, data$fIn.tr, data$sOut.tr, extargs, base,
                                                   ind.vl, data$sIn.vl, data$fIn.vl, data$sOut.vl, active))
        } else {
          cat("\n")
          sub_result[[j]] <- fitNtests_ACO(ants[i,], data$sIn.tr, data$fIn.tr, data$sOut.tr, extargs, base,
                                             ind.vl, data$sIn.vl, data$fIn.vl, data$sOut.vl, active)
        }

        # check if we are still on time
        dt <- difftime(Sys.time(), time.str, units = 'secs')
        if (dt >= time.lim) {
          if (j < n.rep) {
            sub_result[(j+1):n.rep] <- NULL
          }
          break
        }
      }
      # extract complete evaluations
      done <- which(!sapply(sub_result, is.null))
      argsList <- lapply(sub_result[done], `[[`, 1)
      modelList <- lapply(sub_result[done], `[[`, 2)
      fitnessVec <- sapply(sub_result[done], `[[`, 3)

      # identify crashes and usable models
      ids.cr <- which(is.na(fitnessVec))
      ids.ok <- which(!is.na(fitnessVec))

      # remove crashes and compute average model fitness
      if (length(ids.ok) > 0) {
        fitness <- mean(fitnessVec[ids.ok])
        piv <- sample.vec(which(fitnessVec == max(fitnessVec)))
        args <- argsList[[piv]]
        model <- modelList[[piv]]
        result[[i]] <- list(args, model, fitness)
      } else {
        result[[i]] <- NULL
      }
      setTxtProgressBar(pb, i)

      dt <- difftime(Sys.time(), time.str, units = 'secs')
      if (dt >= time.lim) {
        if (i < n.pop) {
          result[(i+1):n.pop] <- NULL
          setTxtProgressBar(pb, n.pop)
        }
        break
      }
    }
    close(pb)
  }
  return(result)
}
# ==========================================================================================================



# ==========================================================================================================
# Function that makes the actual call to fgpm in order to build each model sent by the ants
# ==========================================================================================================
fitNtests_ACO <- function(ant, sIn, fIn, sOut, extargs, base,
                            ind.vl = NULL, sIn.vl = NULL, fIn.vl = NULL, sOut.vl = NULL, active = NULL) {
  b.model <- b.args <- NULL
  args <- formatSol_ACO(ant, sIn, fIn, base)
  poterr <- tryCatch(
    {
      model <- fgpm(sIn = args$sIn, fIn = args$fIn, sOut = sOut, kerType = args$kerType,
                     f_disType = args$f_disType, f_pdims = args$f_pdims, f_basType = args$f_basType,
                     nugget = extargs$nugget, n.starts = extargs$n.starts, n.presample = extargs$n.presample)
    },
    error = function(e) e
  )

  # if the model was succesfully built, compute model fitness
  if (!inherits(poterr, "error")) {
    fitness <- getFitness(model, sIn.vl, fIn.vl, sOut.vl, active)
  } else {
    model <- poterr
    fitness <- NA
  }
  return(list(args[-(1:2)], model, fitness))
}
# ==========================================================================================================



# ==========================================================================================================
# Function to prepare input data structures for prediction based on a fgpm arguments
# ==========================================================================================================
#' @title Preparation of inputs for predictions based on an fgpm modelCall
#' @description This function prepared input data structures according to the active inputs specified by a
#'   \code{"\linkS4class{modelCall}"} object. This function is intended to easily adapt the data structures
#'   to the requirements of a specific model delivered by the model factory function
#'   \link[funGp]{fgpm_factory}.
#'
#' @param sIn.pr sIn.pr an optional matrix of scalar input coordinates at which the output values should be
#'   predicted. Each column is interpreted as a scalar input variable and each row as a coordinate.
#'   Either scalar input coordinates (sIn.pr), functional input coordinates (fIn.pr), or both must be provided.
#'   The \code{"\linkS4class{modelCall}"} object provided through args will lead the extraction of only the
#'   active scalar inputs in the model.
#' @param fIn.pr an optional list of functional input coordinates at which the output values should be
#'   predicted. Each element of the list is interpreted as a functional input variable. Every functional input
#'   variable should be provided as a matrix with one curve per row. Either scalar input coordinates (sIn.pr),
#'   functional input coordinates (fIn.pr), or both must be provided. The \code{"\linkS4class{modelCall}"}
#'   object provided through args will lead the extraction of only the active functional inputs in the model.
#' @param args an object of class \code{"\linkS4class{modelCall}"}, which specifies the set of active
#'   scalar and functional inputs.
#'
#' @return An object of class \code{"list"}, containing the input data structures with only the active inputs
#'   specified by \emph{args}.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @references Betancourt, J., Bachoc, F., and Klein, T. (2020),
#' R Package Manual: "Gaussian Process Regression for Scalar and Functional Inputs with funGp - The in-depth tour".
#' \emph{RISCOPE project}.
#' \href{https://hal.archives-ouvertes.fr/hal-02536624}{[HAL]}
#'
#' @seealso \strong{*} \link[funGp]{predict} for predictions based on a funGp model;
#' @seealso \strong{*} \link[funGp]{fgpm} for creation of a funGp model;
#' @seealso \strong{*} \link[funGp]{fgpm_factory} for funGp heuristic model selection.
#'
#' @importFrom qdapRegex rm_between
#' @export
format4pred <- function(sIn.pr = NULL, fIn.pr = NULL, args) {
  s.on <- f.on <- F
  sIn.4mat <- fIn.4mat <- NULL

  key <- args@string # extract the model call
  if (grepl("sIn", key)) { # is any scalar input active?
    s.on <- T
  }
  if (grepl("fIn", key)) { # is any functional input active?
    f.on <- T
  }

  if (all(s.on, f.on)) { # there both, scalar and funtional inputs active
    # formatting scalar inputs ----------------------------------------------------
    x <- rm_between(key, "sIn = ", ", fIn", extract = TRUE)[[1]] # extract the fragment of the call related to the scalar inputs
    if (grepl("\\d", x)) { # does the fragment include input indices?
      if (grepl("\\(", x)) { # are they grouped in a vector?
        sIn.4mat <- sIn.pr[, as.numeric(unlist(regmatches(x, gregexpr("[[:digit:]]+", x)))), drop=F]
      } else if (grepl("drop", x)) {
        sIn.4mat <- sIn.pr[, as.numeric(unlist(regmatches(x, gregexpr("[[:digit:]]+", x)))), drop=F]
      } else {
        v <- as.numeric(unlist(regmatches(x, gregexpr("[[:digit:]]+", x))))
        sIn.4mat <- sIn.pr[, v[1]:v[2]]
      }

    } else {
      sIn.4mat <- sIn.pr
    }

    # formatting functional inputs ----------------------------------------------------
    x <- rm_between(key, "fIn = ", ", sOut", extract = TRUE)[[1]] # extract the fragment of the call related to the functional inputs
    if (grepl("\\d", x)) { # does the fragment include input indices?
      if (grepl("\\(", x)) { # are they grouped in a vector?
        fIn.4mat <- fIn.pr[[as.numeric(unlist(regmatches(x, gregexpr("[[:digit:]]+", x))))]]
      } else {
        fIn.4mat <- fIn.pr[as.numeric(rm_between(x, "[", "]", extract = TRUE)[[1]])]
      }
    } else {
      fIn.4mat <- fIn.pr
    }

  } else if (s.on) { # there are only scalar inputs active
    # formatting scalar inputs ----------------------------------------------------
    x <- rm_between(key, "sIn = ", ", sOut", extract = TRUE)[[1]] # extract the fragment of the call related to the scalar inputs
    if (grepl("\\d", x)) { # does the fragment include input indices?
      if (grepl("\\(", x)) { # are they grouped in a vector?
        sIn.4mat <- sIn.pr[, as.numeric(unlist(regmatches(x, gregexpr("[[:digit:]]+", x)))), drop=F]
      } else if (grepl("drop", x)) {
        sIn.4mat <- sIn.pr[, as.numeric(unlist(regmatches(x, gregexpr("[[:digit:]]+", x)))), drop=F]
      } else {
        v <- as.numeric(unlist(regmatches(x, gregexpr("[[:digit:]]+", x))))
        sIn.4mat <- sIn.pr[, v[1]:v[2]]
      }

    } else {
      sIn.4mat <- sIn.pr
    }

  } else if (f.on) { # there are only functional inputs active
    # formatting functional inputs ----------------------------------------------------
    x <- rm_between(key, "fIn = ", ", sOut", extract = TRUE)[[1]] # extract the fragment of the call related to the functional inputs
    if (grepl("\\d", x)) { # does the fragment include input indices?
      if (grepl("\\(", x)) { # are they grouped in a vector?
        fIn.4mat <- fIn.pr[[as.numeric(unlist(regmatches(x, gregexpr("[[:digit:]]+", x))))]]
      } else {
        fIn.4mat <- fIn.pr[as.numeric(rm_between(x, "[", "]", extract = TRUE)[[1]])]
      }
    } else {
      fIn.4mat <- fIn.pr
    }

  } else return(NULL)

  return(list(sIn.pr = sIn.4mat, fIn.pr = fIn.4mat))
}
# ==========================================================================================================
