# ==========================================================================================================
# Admin and preprocessing
# ==========================================================================================================

# Method to plot a funGp model
# ----------------------------------------------------------------------------------------------------------
master_ACO <- function(sIn, fIn, sOut, solspace, setup) {
  # set heuristic parameters based on defaults and user specifications
  param <- setParams_ACO(setup)

  # set pheromones and visibility considering the solution space and initial parameters
  env <- setEnvir_ACO(solspace, param)

  # perform exploration
  res <- run_ACO(sIn, fIn, sOut, param, env, solspace$sp.base)

  # correct the howCalled of best model
  res$model@howCalled@string <- getCall_ACO(sIn, fIn, res$args)

  # fill antLog
  res$log <- getLog_ACO(sIn, fIn, res$log.vec, res$log.fitness, solspace$sp.base)

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
# n.lbest: number of best ants to be used for the global pheromone update
# rho.g: global reinforcement coefficient: the larger it is, the more influence best ants have on each update
# =====================================================================================================
setParams_ACO <- function(setup) {
  if (!is.null(setup$n.gen)) n.gen <- setup$n.gen else n.gen <- 10
  if (!is.null(setup$n.pop)) n.pop <- setup$n.pop else n.pop <- 50
  if (!is.null(setup$tao0)) tao0 <- setup$tao0 else tao0 <- 10^-8
  if (!is.null(setup$vis.s)) vis.s <- setup$vis.s else vis.s <- .7
  if (!is.null(setup$vis.f)) vis.f <- setup$vis.f else vis.f <- .7
  if (!is.null(setup$dec.f)) dec.f <- setup$dec.f else dec.f <- .4
  if (!is.null(setup$q0)) q0 <- setup$q0 else q0 <- .95
  if (!is.null(setup$alp)) alp <- setup$alp else alp <- 10
  if (!is.null(setup$bet)) bet <- setup$bet else bet <- 2
  if (!is.null(setup$rho.l)) rho.l <- setup$rho.l else rho.l <- 0
  if (!is.null(setup$dt.l)) dt.l <- setup$dt.l else dt.l <- tao0
  if (!is.null(setup$u.gbest)) u.gbest <- setup$u.gbest else u.gbest <- F
  if (!is.null(setup$n.lbest)) n.lbest <- setup$n.lbest else n.lbest <- 1
  if (!is.null(setup$rho.g)) rho.g <- setup$rho.g else rho.g <- 1

  params <- list(n.gen = n.gen, n.pop = n.pop, tao0 = tao0, vis.s = vis.s, vis.f = vis.f,
                 dec.f = dec.f, q0 = q0, alp = alp, bet = bet, rho.l = rho.l, dt.l = dt.l,
                 u.gbest = u.gbest, n.lbest = n.lbest, rho.g = rho.g)

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

getCall_ACO <- function(sIn, fIn, args) {
  # recover input dimensions
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0

  # identify used scalar inputs
  if (!is.null(args$sIn)) {
    s.used <- which(apply(sIn, 2, function(v) c.vecInMat_Match(args$sIn, v)))
    if (length(s.used) == ncol(sIn)) {
      s.str <- "sIn"
    } else if (length(s.used) == 1) {
      s.str <- paste("sIn[,", s.used, "]", sep = "")
    } else if (all(diff(s.used) == 1)) {
      s.str <- paste("sIn[,", paste(range(s.used), collapse = ":"), "]", sep = "")
    } else {
      s.str <- paste("sIn[,c(", paste(s.used, collapse = ","), ")]", sep = "")
    }
  }

  # identify used functional inputs
  if (!is.null(args$fIn)) {
    f.used <- which(sapply(fIn, function(M) matInList_Match(args$fIn, M)))
    if (length(f.used) == length(fIn)) {
      f.str <- "fIn"
    } else if (length(f.used) == 1) {
      f.str <- paste("fIn[[", f.used, "]]", sep = "")
    } else if (all(diff(f.used) == 1)) {
      f.str <- paste("fIn[[", paste(range(f.used), collapse = ":"), "]]", sep = "")
    } else {
      f.str <- paste("fIn[[c(", paste(f.used, collapse = ","), ")]]", sep = "")
    }

    # prepare f_distype string
    if (length(args$f_disType) == 1) {
      f_disType.str <- paste('"', args$f_disType, '"', sep = "")
    } else {
      f_disType.str <- paste('c("', paste(args$f_disType, collapse = '", "'), '")', sep = "")
    }

    # prepare f_pdims string
    if (length(args$f_pdims) == 1) {
      f_pdims.str <- args$f_pdims
    } else {
      f_pdims.str <- paste('c(', paste(args$f_pdims, collapse = ', '), ')', sep = "")
    }

    # prepare f_basType string
    if (length(args$f_basType) == 1) {
      f_basType.str <- paste('"', args$f_basType, '"', sep = "")
    } else {
      f_basType.str <- paste('c("', paste(args$f_basType, collapse = '", "'), '")', sep = "")
    }
  }

  # prepare kerType string
  kerType.str <- paste('"', args$kerType, '"', sep = "")

  # merge strings to produce model call
  if (all(ds > 0, df > 0)) {
    modcall <- paste("funGp(sIn = ", s.str, ", fIn = ", f.str, ", sOut = sOut, kerType = ", kerType.str,
                     ", f_disType = ", f_disType.str, ", f_pdims = ", f_pdims.str, ", f_basType = ", f_basType.str, ")", sep = "")
  } else if (df > 0) {
    modcall <- paste("funGp(fIn = ", f.str, ", sOut = sOut, kerType = ", kerType.str,
                     ", f_disType = ", f_disType.str, ", f_pdims = ", f_pdims.str, ", f_basType = ", f_basType.str, ")", sep = "")
  } else {
    modcall <- paste("funGp(sIn = ", s.str, ", sOut = sOut, kerType = ", kerType.str, ")", sep = "")
  }

  return(modcall)
}

getLog_ACO <- function(sIn, fIn, log.vec, log.fitness, base) {
  mylog <- new("antsLog")
  args <- list()
  for (i in 1:nrow(log.vec)) {
    mc <- new("modelCall")
    mc@string <- getCall_ACO(sIn, fIn, formatSol_ACO(log.vec[i,], sIn, fIn, base))
    args[[i]] <- mc
  }
  mylog@args <- args
  mylog@sols <- data.frame(log.vec)
  mylog@fitness <- log.fitness
  return(mylog)
}
