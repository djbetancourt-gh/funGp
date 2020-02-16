# ==========================================================================================================
# Admin and preprocessing
# ==========================================================================================================

# Method to plot a funGp model
# ----------------------------------------------------------------------------------------------------------
master_ACO <- function(sIn, fIn, sOut, ind.vl, solspace, setup, extargs, start.time, time.lim, quietly) {
  # set heuristic parameters based on defaults and user specifications
  param <- setParams_ACO(setup)

  # set pheromones and visibility considering the solution space and initial parameters
  env <- setEnvir_ACO(solspace, param)

  # perform exploration
  res <- run_ACO(sIn, fIn, sOut, ind.vl, param, env, solspace$sp.base, extargs, start.time, time.lim, quietly)

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

# =====================================================================================================
# ACO parameters checklist
# =====================================================================================================
# <---> population factors
# n.gen: number of generations to run
# n.pop: number of ants per generation
#
# <---> initial pheromones and visibility
# tao0: initial pheromones (except for cases of already selected factors). should be decimal!
# vis.s: probability of activating any scalar input
# vis.f: probability of activating any functional input
# dec.f: decay rate for visbility of dimension k for the functional inputs (loss: exp(-dec.f * k))
#
# <---> transition rules
# q0: probability of using the greedy rule (more exploitation), otherwise the random proportional rule is used (more exploration)
# alp: influence of the pheromones in random rule. !!!!!!!!! tao^-alp should not give Inf (preferably alp in [0,10])
# bet: influence of the visibility random rule. !!!!!!!!!! nab^-bet should not give Inf (preferably bet in [0,10])
#
# <---> local pheromone update
# rho.l: local evaporation rate: the lower it is, the more time learning is preserved
# dt.l: compensation factor for local update
#
# <---> global pheromone update
# u.gbest: should the global best ant be used for the global update?
# n.gbest: number of best ants to be used for the global pheromone update
# rho.g: global reinforcement coefficient: the larger it is, the more influence best ants have on each update
# =====================================================================================================
setParams_ACO <- function(setup) {
  if (!is.null(setup$n.gen)) n.gen <- setup$n.gen else n.gen <- 5
  if (!is.null(setup$n.pop)) n.pop <- setup$n.pop else n.pop <- 10
  if (!is.null(setup$tao0)) tao0 <- setup$tao0 else tao0 <- 10^-8
  if (!is.null(setup$vis.s)) vis.s <- setup$vis.s else vis.s <- .7
  if (!is.null(setup$vis.f)) vis.f <- setup$vis.f else vis.f <- .7
  if (!is.null(setup$dec.f)) dec.f <- setup$dec.f else dec.f <- .4
  if (!is.null(setup$q0)) q0 <- setup$q0 else q0 <- .9
  if (!is.null(setup$alp)) alp <- setup$alp else alp <- 1
  if (!is.null(setup$bet)) bet <- setup$bet else bet <- 2
  if (!is.null(setup$rho.l)) rho.l <- setup$rho.l else rho.l <- .1
  if (!is.null(setup$dt.l)) dt.l <- setup$dt.l else dt.l <- tao0
  if (!is.null(setup$u.gbest)) u.gbest <- setup$u.gbest else u.gbest <- F
  if (!is.null(setup$n.gbest)) n.gbest <- setup$n.gbest else n.gbest <- 1
  if (!is.null(setup$rho.g)) rho.g <- setup$rho.g else rho.g <- .1

  params <- list(n.gen = n.gen, n.pop = n.pop, tao0 = tao0, vis.s = vis.s, vis.f = vis.f,
                 dec.f = dec.f, q0 = q0, alp = alp, bet = bet, rho.l = rho.l, dt.l = dt.l,
                 u.gbest = u.gbest, n.gbest = n.gbest, rho.g = rho.g)

  return(params)
}

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
  vis.s <- param$vis.s
  vis.f <- param$vis.f
  dec.f <- param$dec.f

  # extract input dimensions
  ds <- sp.user$ds
  df <- sp.user$df

  # create visibility and pheromones lists
  visib <- phero <- list()
  layers <- c()

  if (ds > 0) {
    # set up visibility and pheromones related to scalar inputs
    if (s.state.u[1] == 0) {
      phero[[1]] <- matrix(tao0, ncol = 2) # set up initial pheromones
      visib[[1]] <- matrix(c(vis.s, (1 - vis.s)), ncol = 2) # set up visibility
    } else {
      phero[[1]] <- matrix(c(1, 0), ncol = 2) # set up initial pheromones
      visib[[1]] <- matrix(c(1, 0), ncol = 2) # set up visibility
    }
    rownames(visib[[1]]) <- rownames(phero[[1]]) <- "Orig"
    colnames(visib[[1]]) <- colnames(phero[[1]]) <- c("Active", "Inactive")
    layers <- c(layers, paste("State X", 1, sep =""))

    if (ds > 1) {
      for (i in 2:ds){
        # state of scalar input i
        # ____________________________________________________________________________________
        if (s.state.u[i] == 0) {
          phero[[i]] <- matrix(tao0, nrow = 2, ncol = 2) # set up initial pheromones
          visib[[i]] <- matrix(rep(c(vis.s, (1 - vis.s)), 2), nrow = 2, byrow = T) # set up visibility
        } else {
          phero[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T) # set up initial pheromones
          visib[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T) # set up visibility
        }
        # set up colnames
        colnames(phero[[i]]) <- rownames(phero[[i]]) <- c("Active", "Inactive")
        colnames(visib[[i]]) <- rownames(visib[[i]]) <- c("Active", "Inactive")
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
      i <- length(visib) + 1
      if (all(i == 1, j == 1)) {
        if (f.state.u[j] == 0) { # free, distribute pheromones and assign default visibilities
          phero[[i]] <- matrix(tao0, ncol = 2) # set up initial pheromones
          visib[[i]] <- matrix(c(vis.f, (1 - vis.f)), ncol = 2) # set up visibility

        } else { # fixed, put all load in first colum
          phero[[i]] <- matrix(c(1, 0), ncol = 2) # set up initial pheromones
          visib[[i]] <- matrix(c(1, 0), ncol = 2) # set up visibility
        }
        rownames(visib[[i]]) <- rownames(phero[[i]]) <- "Orig"
        colnames(visib[[i]]) <- colnames(phero[[i]]) <- c("Active", "Inactive")

      } else {
        if (f.state.u[j] == 0) { # free, distribute pheromones and assign default visibilities
          phero[[i]] <- matrix(tao0, nrow = 2, ncol = 2)
          visib[[i]] <- matrix(rep(c(vis.f, (1-vis.f)), 2), nrow = 2, byrow = T)

        } else { # fixed, put all load in first colum
          phero[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
          visib[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
        }
        rownames(visib[[i]]) <- rownames(phero[[i]]) <- colnames(phero[[i-1]])
        colnames(visib[[i]]) <- colnames(phero[[i]]) <- c("Active", "Inactive")
      }
      layers <- c(layers, paste("State F", j, sep =""))
      # ____________________________________________________________________________________


      # distance for functional input j
      # ____________________________________________________________________________________
      i <- length(visib) + 1
      nr <- ncol(visib[[i-1]])
      # <--> identify active levels
      v1 <- v2 <- rep(0, length(f.dist.0[[j]]))
      t <- unname(sapply(f.dist.u[[j]], function(x) which(f.dist.0[[j]] == x)))
      # <--> assign initial pheromones
      if (length(t) == 1)  v1[t] <- 1 else v1[t] <- tao0
      phero[[i]] <- matrix(rep(v1, nr), nrow = nr, byrow = T)
      # <--> distribute visbility evenly
      v2[t] <- 1
      visib[[i]] <- matrix(rep(v2/sum(v2), nr), nrow = nr, byrow = T)
      # <--> names
      rownames(visib[[i]]) <- rownames(phero[[i]]) <- colnames(phero[[i-1]])
      colnames(visib[[i]]) <- colnames(phero[[i]]) <- f.dist.0[[j]]
      layers <- c(layers, paste("Dist. FI", j, sep =""))
      # ____________________________________________________________________________________


      # dimension for functional input j
      # ____________________________________________________________________________________
      i <- length(visib) + 1
      nr <- ncol(visib[[i-1]])
      # <--> identify active levels
      v <- vgro <- vind <- rep(0, length(f.dims.0[[j]])) # zero vector with as many elements as potential dimensions in the original set
      t <- unname(sapply(f.dims.u[[j]], function(x) which(f.dims.0[[j]] == x))) # locations of dimensions still active
      # <--> assign initial pheromones
      if (length(t) == 1)  v[t] <- 1 else v[t] <- tao0
      phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
      # <--> visbility vith L2_bygroup distance (distribute evenly)
      if (length(t) == 1)  vgro[t] <- 1 else vgro[t] <- 1/length(t)
      # <--> visbility vith L2_byindex distance (distribute according to loss function)
      if (0 %in% f.dims.u[[j]]) l <- exp(-dec.f * c(max(f.dims.u[[j]]), f.dims.u[[j]][-1])) else l <- exp(-dec.f * f.dims.u[[j]])
      if (length(t) == 1)  vind[t] <- 1 else vind[t] <- l
      # <--> wrap visibility
      visib[[i]] <- rbind(vgro, vind)
      # <--> names
      # <--> names
      rownames(visib[[i]]) <- rownames(phero[[i]]) <- colnames(phero[[i-1]])
      colnames(visib[[i]]) <- colnames(phero[[i]]) <- f.dims.0[[j]]
      layers <- c(layers, paste("Dim. FI", j, sep =""))
      # ____________________________________________________________________________________


      # basis for functional input j
      # ____________________________________________________________________________________
      i <- length(visib) + 1
      nr <- ncol(visib[[i-1]])
      # <--> identify active levels
      v1 <- v2 <- rep(0, length(f.bas.0[[j]])) # zero vector with as many elements as basis families in the original set
      t <- unname(sapply(f.bas.u[[j]], function(x) which(f.bas.0[[j]] == x))) # locations of basis still active
      # <--> assign initial pheromones
      if (length(t) == 1)  v1[t] <- 1 else v1[t] <- tao0
      phero[[i]] <- matrix(rep(v1, nr), nrow = nr, byrow = T)
      # <--> distribute visbility evenly
      v2[t] <- 1
      visib[[i]] <- matrix(rep(v2/sum(v2), nr), nrow = nr, byrow = T)
      # <--> names
      colnames(visib[[i]]) <- colnames(phero[[i]]) <- f.bas.0[[j]]
      rownames(visib[[i]]) <- rownames(phero[[i]]) <- f.dims.0[[j]]
      layers <- c(layers, paste("Basis FI", j, sep =""))
      # ____________________________________________________________________________________
    }
  }

  # kernel function
  # ____________________________________________________________________________________
  i <- length(visib) + 1
  nr <- ncol(visib[[i-1]])
  # <--> identify active levels
  v1 <- v2 <- rep(0, length(k.type.0)) # zero vector with as many elements as kernel functions in the original set
  t <- unname(sapply(k.type.u, function(x) which(k.type.0 == x))) # locations of kernels still active
  # <--> assign initial pheromones
  if (length(t) == 1)  v1[t] <- 1 else v1[t] <- tao0
  phero[[i]] <- matrix(rep(v1, nr), nrow = nr, byrow = T)
  # <--> distribute visbility evenly
  v2[t] <- 1
  visib[[i]] <- matrix(rep(v2/sum(v2), nr), nrow = nr, byrow = T)
  # <--> names
  rownames(visib[[i]]) <- rownames(phero[[i]]) <- colnames(phero[[i-1]])
  colnames(visib[[i]]) <- colnames(phero[[i]]) <- k.type.0
  layers <- c(layers, "Kernel")
  # ____________________________________________________________________________________

  # set names to visibility and pheromones lists
  names(visib) <- names(phero) <- layers

  return(list(phero = phero, visib = visib))
}

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
      # print(paste("s.active:", s.active.ac)) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Remove this line
      sIn.ac <- sIn[,s.active.ac, drop = F]
    } else {
      sIn.ac <- NULL
    }
  } else {
    sIn.ac <- NULL
  }
  piv <- ds

  defargs <- formals(funGp) # in case they are required below
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
        # browser()
        f.dist.ac[[j]] <- defargs$f_disType
        f.dims.ac[[j]] <- defargs$f_pdims
        f.bas.ac[[j]] <- defargs$f_basType
        piv <- piv + 2
      } else {
        f.dist.ac[[j]] <- f.dist[[j]][ant[piv]]
        piv <- piv + 1
        # if (ant[piv] == 1) f.dims.ac[[j]] <- 0 else f.dims.ac[[j]] <- ant[piv] - 1
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

getActiveIn_ACO <- function(ant, sIn, fIn, base) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  # recover base components
  s.state <- base$s.state
  f.state <- base$f.state

  # index of active scalar inputs
  if (ds > 0) s.active.ac <- unname(which(ant[1:ds] == 1)) else s.active.ac <- NULL

  # index of active functional inputs
  if (df > 0) f.active.ac <- unname(which(ant[grepl("State F", names(ant))] == 1)) else f.active.ac <- NULL

  return(list(s.active = s.active.ac, f.active = f.active.ac))
}

getCall_ACO <- function(ant, sIn, fIn, args, base, extargs) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  # identify used inputs
  in.used <- getActiveIn_ACO(ant, sIn, fIn, base)

  # set string for scalar inputs
  if (!is.null(args$sIn)) {
    # s.used <- which(apply(sIn, 2, function(v) c.vecInMat_Match(args$sIn, v)))
    s.used <- in.used$s.active
    if (length(s.used) == ncol(sIn)) {
      s.str <- "sIn = sIn"
    } else if (length(s.used) == 1) {
      s.str <- paste("sIn = sIn[,", s.used, "]", sep = "")
    } else if (all(diff(s.used) == 1)) {
      s.str <- paste("sIn = sIn[,", paste(range(s.used), collapse = ":"), "]", sep = "")
    } else {
      s.str <- paste("sIn = sIn[,c(", paste(s.used, collapse = ","), ")]", sep = "")
    }
  } else {
    s.str <- NULL
  }

  # set string for functional inputs
  if (!is.null(args$fIn)) {
    # f.used <- which(sapply(fIn, function(M) matInList_Match(args$fIn, M)))
    f.used <- in.used$f.active
    if (length(f.used) == length(fIn)) {
      f.str <- "fIn = fIn"
    } else if (length(f.used) == 1) {
      f.str <- paste("fIn = fIn[[", f.used, "]]", sep = "")
    } else if (all(diff(f.used) == 1)) {
      f.str <- paste("fIn = fIn[[", paste(range(f.used), collapse = ":"), "]]", sep = "")
    } else {
      f.str <- paste("fIn = fIn[[c(", paste(f.used, collapse = ","), ")]]", sep = "")
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

  # prepare kerType string
  kerType.str <- paste('kerType = "', args$kerType, '"', sep = "")

  # set strings for extra arguments: nugget, n.starts, n.presample
  defargs <- formals(funGp)
  if (extargs$nugget != defargs$nugget) nugg.str <- paste("nugget = ", extargs$nugget, sep = "") else nugg.str <- NULL
  if (extargs$n.starts != defargs$n.starts) n.starts.str <- paste("n.starts = ", extargs$n.starts, sep = "") else n.starts.str <- NULL
  if (extargs$n.presample != defargs$n.presample) {
    n.presample.str <- paste("n.presample = ", extargs$n.presample, sep = "")
  } else {
    n.presample.str <- NULL
  }

  # merge strings to produce model call
  full.str <- c(s.str, f.str, f_disType.str, f_pdims.str, f_basType.str, kerType.str, nugg.str, n.starts.str, n.presample.str)
  modcall <- paste("funGp(", paste(full.str, collapse = ", "), ")", sep = "")

  return(modcall)
}

vec2DFrame_ACO <- function(sol.vec, sIn, fIn) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  names <- vals <- rep(0,length(sol.vec))
  if (ds > 0) {
    for (i in 1:ds) {
      names[i] <- paste("State_X", i, sep = "")
      # vals[i] <- c("Active", "Inactive")[sol.vec[i]]
      vals[i] <- c("On", "Off")[sol.vec[i]]
    }
    piv <- ds
  }

  if (df > 0) {
    piv <- max(0, ds)
    for (i in 1:df) {
      piv <- piv + 1
      names[piv] <- paste("State_F", i, sep = "")
      # vals[piv] <- c("Active", "Inactive")[sol.vec[piv]]
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
