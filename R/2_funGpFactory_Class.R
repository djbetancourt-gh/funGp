# ==========================================================================================================
# Class for factories to produce structre-optimized funGp models
# ==========================================================================================================



# ==========================================================================================================
# Developer oriented methods
# ==========================================================================================================

# Constructor of the class
# ----------------------------------------------------------------------------------------------------------
#' @title Class: Fill!!!!
#' @description Fill this!!!!!!!!!
#'
#' @slot ppp Object of class \code{"character"}. Fill!!!!!!!!!!
#'
#' @rdname factory-class
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setClass("funGpFactory",
         representation(
           ppp = "character"              # a ppp
         ),
         validity = function(object) {T})
# ----------------------------------------------------------------------------------------------------------



# funGp master function: used for construction and training of a funGP model
# ----------------------------------------------------------------------------------------------------------
#' @title Also Fill!!
#' @description Fill!!.
#' @param sIn a matrix of scalar input values to train the model. Each column must match an input variable
#' @param fIn a list of functional inputs values to fit the model. Each element of the list must contain a
#' matrix.
#' @param sOut a vector (or 1-column matrix) containing the values of the scalar output at the training
#' points.
#' @param tr.ind matrix of training indices.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
symbolicfunGp <- function(sIn = NULL, fIn = NULL, sOut, tr.ind = NULL) {
# symbolicfunGp <- function(sIn = NULL, fIn = NULL, sOut, tr.ind = NULL, q0, alp, bet, rho.l, rho.g) {
  print("Hey, lets optimize that modal!")

  if (!is.null(sIn)) {sIn <- as.matrix(sIn); ds <- ncol(sIn)} else ds <- 0
  if (!is.null(fIn)) {
    if (!is.list(fIn)) {
      fIn <- list(fIn)
    }
    df <- length(fIn)
  } else {
    df <- 0
  }

  # define template with base structures and values for a full experiment
  s.state <- rep(0, ds) # 0 = free, 1 = fixed
  f.state <- rep(0, df) # 0 = free, 1 = fixed
  f.dims <- lapply(sapply(fIn, ncol), function(k) 0:k) # sequence of potential dimensions p for each variable. 0 means no proj
  f.dist <- rep(list(c("L2_bygroup", "L2_byindex")), df) # all available distances at the time
  f.fam <- rep(list(c("B-splines", "PCA")), df) # all available basis families for functions at the time
  k.type <- c("gaussian", "matern5_2", "matern3_2")
  base <- list(s.state = s.state, f.state = f.state, f.dims = f.dims, f.dist = f.dist, f.fam = f.fam, k.type = k.type)

  # make a copy for the user to modify
  copy <- base

  # user inputs, remove this at the end !!!!!!!!!!!!!!!!!!!!!!!!!
  s_fixed <- c(1) # keep X1 always active
  f_fixed <- NULL # Do not keep any functional input always active, let them free
  f_fixdim <- matrix(c(2,4), ncol = 1) # set dimension of F2 at 4
  f_maxdim <- matrix(c(1,5), ncol = 1) # set max dimension of F1 at 5
  f_setdist <- list("2" = c("L2_byindex")) # only test L2_index distance for F2
  f_setfam <- list("1" = c("B-splines")) # only B-splines projection for F1
  k_settype <- c("matern5_2", "matern3_2") # test only matern kernels

  # update copies based on user inputs into
  clean <- updateSpace(base = base, copy = copy, ds, df,
                       s_fixed = s_fixed, f_fixed = f_fixed, f_fixdim = f_fixdim, f_maxdim = f_maxdim,
                       f_setdist = f_setdist, f_setfam = f_setfam, k_settype = k_settype)

  # print template and actual exeriment
  # checkAndPrint(ds, df, base)
  # checkAndPrint(ds, df, clean)

  # set up heuristic parameters
  #___________________________________________________________________________________________
  # <---> population factors
  n.gen <- 10 # number of generations to run
  n.pop <- 5 # 30 # number of ants per generation

  # <---> initial pheromones and visibility
  tao0 <- 10^-8 # initial pheromones (except for cases of already selected factors). should be decimal!
  vis.s <- .7 # probability of activating any scalar input
  vis.f <- .7 # probability of activating any functional input
  dec.f <- .4 # decay rate for visbility of dimension k for the functional inputs (loss: exp(-dec.f * k))

  # <---> transition rules
  q0 <- 1 #.7 # probability of using the greedy rule (more exploitation), otherwise the random proportional rule is used (more exploration)
  alp <- 1 # influence of the pheromones in random rule. !!!!!!!!! tao^-alp should not give Inf (preferably alp in [0,10])
  bet <- 2 # influence of the visibility random rule. !!!!!!!!!! nab^-bet should not give Inf (preferably bet in [0,10])

  # <---> local pheromone update
  rho.l <- .1 # local evaporation rate: the lower it is, the more time learning is preserved
  dt.l <- tao0 # compensation factor for local update

  # <---> global pheromone update
  u.gbest <- F # should the global best ant be used for the global update?
  n.lbest <- 1 # number of best ants to be used for the global pheromone update
  rho.g <- .1 # global reinforcement coefficient: the larger it is, the more influence best ants have on each update
  #___________________________________________________________________________________________

  # set up initial environment conditions (visibility and pheromones)
  #___________________________________________________________________________________________
  tmp <- setPherAndVis(ds, df, base, clean, tao0, vis.s, vis.f, dec.f)
  phero <- tmp$phero # pheromones
  visib <- tmp$visib # visibility
  #___________________________________________________________________________________________

  # set up heuristic controllers and statistics
  #___________________________________________________________________________________________
  # <---> required controllers
  n.layers <- length(phero) # number of layers in the network/factors in the experiment (some may be fixed constant)

  # <---> for statistics
  evol.fitness <- rep(0, n.gen) # fitness of best ant of each colony
  all.fitness <- matrix(nrow = n.pop, ncol = n.gen)

  # <---> best solution and its fitness
  b.ant <- "just use the mean" # initialize best ant
  b.fitness <- 0 # initialize best fitness
  #___________________________________________________________________________________________

  # browser()
  # run the colony
  for (c.gen in 1:n.gen) {
    start_time <- Sys.time()
    cat(paste("Dispatching colony", c.gen, "\n"))
    pb <- txtProgressBar(min = 0, max = n.pop, style = 3)

    # create a new colony and mark as incomplete
    ants <- matrix(nrow = n.pop, ncol = n.layers)
    colnames(ants) <- names(phero)
    antsPartial <- rep(T, n.pop)
    antsLayers <- rep(0, n.pop)
    antsLevels <- rep(1, n.pop) # everbody starts at the origin (in R is position 1)
    colPartial <- T

    # iterate until the colony is complete
    while (colPartial) {
      # identify incomplete ants
      id.partial <- which(antsPartial)

      # randomly pick an incomplete ant
      id.ant <- sample.vec(id.partial, 1)
      myant <- list(id = id.ant, sol = ants[id.ant,], layer = antsLayers[id.ant], level = antsLevels[id.ant])

      # cat(paste("\nid.ant:", id.ant, "\n"))

      # add step to selected ant based on decision rules
      # print("----------------< Next step")
      if (runif(1) <= q0) {
        # print("greedy")
        antup <- step_greedy(myant, phero, visib, alp, bet, c.gen)
      } else {
        # print("propor")
        antup <- step_propor(myant, phero, visib, alp, bet, c.gen)
      }

      # perform local pheromone update only if the factor was not already fixed
      # print("----------------< Local update")
      if (antup$layer - myant$layer == 1) {
        phero <- local_pheroUpd(phero, myant, antup, rho.l, dt.l)
      }

      # save the ant
      ants[antup$id,] <- antup$sol
      antsLayers[antup$id] <- antup$layer
      antsLevels[antup$id] <- antup$level

      # check if the ant is done
      if (antup$layer == n.layers) {
        antsPartial[antup$id] <- F
        # setTxtProgressBar(pb, sum(!antsPartial))
      }

      # check if the colony is done
      if (all(!antsPartial)) colPartial <- F
    }

    # compute fitness of each ant
    fitness <- rep(0, n.pop)
    for (i in 1:n.pop) {
      # translate ant data into funGp arguments format
      args <- getArgs(ants[i,], ds, df, sIn, fIn, base)

      if (all(is.null(args$sIn), is.null(args$fIn))) browser()

      # build the model
      model <- quiet(funGp(sIn = args$sIn, fIn = args$fIn, sOut = sOut, kerType = args$kerType,
                     f_disType = args$f_disType, f_pdims = args$f_pdims, f_family = args$f_family))

      # compute model fitness #!!!!!!!!!!!!!!!!!! extend to other metrics
      # getFitness(model, sIn.vl, fIn.vl, sOut.vl)
      fitness[i] <- max(getFitness(model),0)
      # print(ants[i,])
      # print(fitness[i])
      setTxtProgressBar(pb, i)
    }
    close(pb)

    # extract ants and fitness for global update
    res <- derssAnts(fitness, n.lbest, ants, u.gbest, c.gen, b.ant, b.fitness)

    # perform global pheromone update
    # print("----------------< Global update")
    phero <- global_pheroUpd(res$ants.up, res$fitness.up, phero, rho.g)

    # save best ant
    if (fitness[res$b.ind[1]] > b.fitness) {
      b.ant <- ants[res$b.ind[1],]
      b.fitness <- fitness[res$b.ind[1]]
    }
    # print(b.ant)

    print(Sys.time() - start_time)
    cat("\n")

    evol.fitness[c.gen] <- b.fitness
    all.fitness[,c.gen] <- fitness

    # plot current best model
    # b.plot <- function(b.ant, b.fitness) {
    #   print(b.fitness)
    #   b.args <- getArgs(b.ant, ds, df, sIn, fIn, base)
    #   model <- funGp(sIn = b.args$sIn, fIn = b.args$fIn, sOut = sOut, kerType = b.args$kerType,
    #                  f_disType = b.args$f_disType, f_pdims = b.args$f_pdims, f_family = b.args$f_family)
    #   plotLOO(model)
    # }

    # if (!is.character(b.ant)) {
    #   b.plot(b.ant, b.fitness)
    # }
  }

  plot(1, type = "n", xlab = "Colony", ylab = "Fitness", xlim = c(1, (n.gen + .3)), ylim = c(0, 1), xaxt = "n")
  axis(1, 1:n.gen)
  for (i in 1:n.gen) {
    points(rep(i, n.pop), all.fitness[,i], pch = 21, bg = scales::alpha("red", .4), col = scales::alpha("red", .4))
    points(i, median(all.fitness[,i]), pch = 21, bg = "blue", col = scales::alpha("blue", .4))
    # legend(x = (i + .2), y = (evol.fitness[i] - .05), legend = format(evol.fitness[i], digits = 2, nsmall = 3), cex = 1,
    #        xjust = 0.5,      # 0.5 means center adjusted
    #        yjust = 0.5,      # 0.5 means center adjusted
    #        x.intersp = -0.5, # adjust character interspacing as you like to effect box width
    #        y.intersp = 0.1,  # adjust character interspacing to effect box height
    #        adj = c(0, 0.5))
  }
  # lines(evol.fitness, lty = 2, col = "blue")
  points(evol.fitness, pch = 21, bg = scales::alpha("blue", 0), col = "black")

  cat("\nAnts are done ;)")
  # return(b.fitness)
}

local_pheroUpd <- function(phero, myant, antup, rho.l, dt.l) {
  # browser()
  o <- myant$level
  d <- antup$level
  if (phero[[antup$layer]][o,d] < 1) {
    # print(paste("phero0:", phero[[antup$layer]][o,d]))
    phero[[antup$layer]][o,d] <- (1 - rho.l) * phero[[antup$layer]][o,d] + rho.l * (dt.l)
    # print(paste("phero1:", phero[[antup$layer]][o,d]))
  } else {
    # print(paste("fixed:", phero[[antup$layer]][o,d]))
  }
  return(phero)
}

global_pheroUpd <- function(b.ants, b.fitness, phero, rho.g) {
  # browser()
  n.best <- nrow(b.ants)
  for (i in 1:n.best) {
    c.ant <- b.ants[i,] # current ant (dynamic during the outer loop)
    c.lev <- 1 # current level (dynamic during the inner loop)
    for (c.lay in 1:length(phero)) { # current layer (dynamic during the inner loop)
      o <- c.lev
      if (grepl("Dim", names(phero)[c.lay])) {
        d <- c.ant[c.lay] + 1 # to correct projection dimension
      } else {
        d <- c.ant[c.lay]
      }

      if (phero[[c.lay]][o,d] < 1) {
        if (all(grepl("State F", names(phero)[c.lay]), c.lev == 2)) {
          # identify the next layer of interest (just before next selection)
          c.lay <- (which(!grepl("FI", names(phero))) > c.lay)[1]
          c.lev <- 1 # current level update
        } else {
          # print(paste("phero0:", phero[[c.lay]][o,d]))
          phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * b.fitness[i] # pheromone update
          c.lev <- d # current level update
          # print(paste("phero1:", phero[[c.lay]][o,d]))
        }
      } else {
        # print(paste("fixed:", phero[[c.lay]][o,d]))
        c.lev <- d
      }
      # c.tao <- phero[[c.lay]][o,d] # current pheromone value
      # p.tao <- (1 - rho.g) * c.tao + rho.g * b.fitness[i] # potential tao # p.tao <- (1 - rho.g) * c.tao + (1 - (1 - rho.g) * c.tao) * b.fitness[i] # potential tao
      # phero[[c.lay]][o,d] <- max(c.tao, p.tao) # pheromone update
    }
  }
  return(phero)
}

step_greedy <- function(myant, phero, visib, alp, bet, c.gen) {
  # browser()
  # copy my ant to make updates on it
  antup <- myant

  # select next step
  layer <- myant$layer + 1
  level <- myant$level
  tao <- unname(phero[[layer]][level,]) # level of pheromones of neighbors
  nab <- unname(visib[[layer]][level,]) # visibility  of neighbors
  # attr <- tao^alp + nab^bet # attractiveness of each neighbor node
  # attr <- tao^(-alp) + nab^(-bet) # attractiveness of each neighbor node
  attr <- alp * tao + bet * nab # attractiveness of each neighbor node
  sel.level <- sample.vec(which(attr == max(attr)), 1)

  # print(paste("sel.level:", sel.level))
  # print(paste("tao:", tao))
  # print(paste("nab:", nab))
  # print(paste("attr:", attr))


  # update ant after selection
  if (grepl("Dim", names(phero)[layer])) {
    antup$sol[layer] <- sel.level - 1 # to correct projection dimension
  } else {
    antup$sol[layer] <- sel.level
  }
  antup$layer <- layer
  antup$level <- sel.level

  # if the selection was to turn off a functional input, then move foward the ant until next input
  if (all(grepl("State F", names(phero)[layer]), sel.level == 2)) {
    browser()
    # identify layers not related to features of a functional input (apart from its state)
    nff.lay <- which(!grepl("FI", names(phero)))
    # identify the next layer of interest (just before next selection)
    tar.lay <- (nff.lay[nff.lay > layer])[1] - 1
    x.skip <- (layer+1):tar.lay # get index of the layers to skip
    antup$sol[x.skip] <- 1 # fill solution at skipped layers (arbitrarily 1 but has no effect)
    antup$layer <- tar.lay # update layer according to skip
    antup$level <- 1 # locate the ant arbitrarily at level 1 (irrelevant for next decision)
  }

  # return the updated ant
  return(antup)
}

step_propor <- function(myant, phero, visib, alp, bet, c.gen) {
  # browser()
  # copy my ant to make updates on it
  antup <- myant

  # select next step
  layer <- myant$layer + 1
  level <- myant$level
  tao <- unname(phero[[layer]][level,]) # level of pheromones of neighbors
  nab <- unname(visib[[layer]][level,]) # visibility  of neighbors
  # attr <- tao^alp + nab^bet # attractiveness of each neighbor node
  attr <- alp * tao + bet * nab # attractiveness of each neighbor node
  sel.level <- sample(1:length(tao), 1, prob = attr/sum(attr))

  # print(paste("sel.level:", sel.level))
  # print(paste("tao:", tao))
  # print(paste("nab:", nab))
  # print(paste("attr:", attr))

  # update ant after selection
  if (grepl("Dim", names(phero)[layer])) {
    antup$sol[layer] <- sel.level - 1 # to correct projection dimension
  } else {
    antup$sol[layer] <- sel.level
  }
  antup$layer <- layer
  antup$level <- sel.level

  # if the selection was to turn off a functional input, then move foward the ant until next input
  if (all(grepl("State F", names(phero)[layer]), sel.level == 2)) {
    # identify layers not related to features of a functional input (apart from its state)
    nff.lay <- which(!grepl("FI", names(phero)))
    # identify the next layer of interest (just before next selection)
    tar.lay <- (nff.lay[nff.lay > layer])[1] - 1
    x.skip <- (layer+1):tar.lay # get index of the layers to skip
    antup$sol[x.skip] <- 1 # fill solution at skipped layers (arbitrarily 1 but has no effect)
    antup$layer <- tar.lay # update layer according to skip
    antup$level <- 1 # locate the ant arbitrarily at level 1 (irrelevant for next decision)
  }

  # return the updated ant
  return(antup)
}


derssAnts <- function(fitness, n.lbest, ants, u.gbest, c.gen, b.ant, b.fitness){
  # browser()
  # identify best n.lbest ants
  b.ind <- order(fitness, decreasing = T)[1:n.lbest]

  # remove duplicates if there is any
  if (n.lbest > 1) {
    u.ind <- b.ind[!duplicated(ants[b.ind,])] # unique best ants
  } else {
    u.ind <- b.ind
  }

  # group best ants and their fitness
  ants.up <- ants[u.ind,,drop = F]
  fitness.up <- fitness[u.ind]

  # add the global best only if required
  if (all(u.gbest, c.gen > 1)) {
    ants.up <- rbind(ants.up, b.ant)
    fitness.up <- c(fitness.up, b.fitness)
  }

  return(list(ants.up = ants.up, fitness.up = fitness.up, b.ind = b.ind))
}


setPherAndVis <- function(ds, df, base, clean, tao0, vis.s, vis.f, dec.f) {
  # recover base components
  s.state.0 <- base$s.state
  f.state.0 <- base$f.state
  f.dims.0 <- base$f.dims
  f.dist.0 <- base$f.dist
  f.fam.0 <- base$f.fam
  k.type.0 <- base$k.type

  # recover clean components
  s.state.c <- clean$s.state
  f.state.c <- clean$f.state
  f.dims.c <- clean$f.dims
  f.dist.c <- clean$f.dist
  f.fam.c <- clean$f.fam
  k.type.c <- clean$k.type

  # create visibility and pheromones lists
  visib <- phero <- list()
  layers <- c()

  # set up initial visibilities and pheromones for scalar inputs
  if (s.state.c[1] == 0) { # free, distribute pheromones and assign default visibilities
    phero[[1]] <- matrix(c(tao0, (1-tao0)), nrow = 1, byrow = T)
    visib[[1]] <- matrix(c(vis.s, (1-vis.s)), nrow = 1, byrow = T)

  } else { # fixed, put all load in active state
    phero[[1]] <- matrix(c(1, 0), nrow = 1, byrow = T)
    visib[[1]] <- matrix(c(1, 0), nrow = 1, byrow = T)
  }
  rownames(visib[[1]]) <- rownames(phero[[1]]) <- "Orig"
  colnames(visib[[1]]) <- colnames(phero[[1]]) <- c("Active", "Inactive")
  layers <- c(layers, paste("State X", 1, sep =""))

  if (ds > 1) {
    for (i in 2:ds) {
      if (s.state.c[i] == 0) { # free, distribute pheromones and assign default visibilities
        phero[[i]] <- matrix(tao0, nrow = 2, ncol = 2) # matrix(rep(c(tao0, (1-tao0)), 2), nrow = 2, byrow = T)
        visib[[i]] <- matrix(rep(c(vis.s, (1-vis.s)), 2), nrow = 2, byrow = T)

      } else { # fixed, put all load in active state
        phero[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
        visib[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
      }
      rownames(visib[[i]]) <- rownames(phero[[i]]) <- c("Active", "Inactive")
      colnames(visib[[i]]) <- colnames(phero[[i]]) <- c("Active", "Inactive")
      layers <- c(layers, paste("State X", i, sep =""))
    }
  }

  # set up visibility and pheromones related to functional inputs
  for (j in 1:df) {
    # state of functional input j
    # ____________________________________________________________________________________
    i <- length(visib) + 1
    if (f.state.c[j] == 0) { # free, distribute pheromones and assign default visibilities
      phero[[i]] <- matrix(tao0, nrow = 2, ncol = 2)
      visib[[i]] <- matrix(rep(c(vis.f, (1-vis.f)), 2), nrow = 2, byrow = T)

    } else { # fixed, put all load in first colum
      phero[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
      visib[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
    }
    rownames(visib[[i]]) <- rownames(phero[[i]]) <- c("Active", "Inactive")
    colnames(visib[[i]]) <- colnames(phero[[i]]) <- c("Active", "Inactive")
    layers <- c(layers, paste("State F", j, sep =""))
    # ____________________________________________________________________________________


    # distance for functional input j
    # ____________________________________________________________________________________
    i <- length(visib) + 1
    nr <- ncol(visib[[i-1]])
    # <--> identify active levels
    v1 <- v2 <- rep(0, length(f.dist.0[[j]]))
    t <- unname(sapply(f.dist.c[[j]], function(x) which(f.dist.0[[j]] == x)))
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
    t <- unname(sapply(f.dims.c[[j]], function(x) which(f.dims.0[[j]] == x))) # locations of dimensions still active
    # <--> assign initial pheromones
    if (length(t) == 1)  v[t] <- 1 else v[t] <- tao0
    phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
    # <--> visbility vith L2_bygroup distance (distribute evenly)
    vgro[t] <- 1/length(t)
    # <--> visbility vith L2_byindex distance (distribute according to loss function)
    if (0 %in% f.dims.c[[j]]) l <- exp(-dec.f * c(max(f.dims.c[[j]]), f.dims.c[[j]][-1])) else l <- exp(-dec.f * f.dims.c[[j]])
    vind[t] <- l
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
    v1 <- v2 <- rep(0, length(f.fam.0[[j]])) # zero vector with as many elements as basis families in the original set
    t <- unname(sapply(f.fam.c[[j]], function(x) which(f.fam.0[[j]] == x))) # locations of basis still active
    # <--> assign initial pheromones
    if (length(t) == 1)  v1[t] <- 1 else v1[t] <- tao0
    phero[[i]] <- matrix(rep(v1, nr), nrow = nr, byrow = T)
    # <--> distribute visbility evenly
    v2[t] <- 1
    visib[[i]] <- matrix(rep(v2/sum(v2), nr), nrow = nr, byrow = T)
    # <--> names
    colnames(visib[[i]]) <- colnames(phero[[i]]) <- f.fam.0[[j]]
    rownames(visib[[i]]) <- rownames(phero[[i]]) <- f.dims.0[[j]]
    layers <- c(layers, paste("Basis FI", j, sep =""))
    # ____________________________________________________________________________________
  }

  # # kernel function
  # ____________________________________________________________________________________
  i <- length(visib) + 1
  nr <- ncol(visib[[i-1]])
  # <--> identify active levels
  v1 <- v2 <- rep(0, length(k.type.0)) # zero vector with as many elements as kernel functions in the original set
  t <- unname(sapply(k.type.c, function(x) which(k.type.0 == x))) # locations of kernels still active
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

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

getOut_loocv <- function(model){
  y_obs <- model@sOut
  R <- tcrossprod(model@preMats$L)/model@kern@varHyp + diag(model@nugget, nrow = model@n.tr, ncol = model@n.tr)
  Rinv <- solve(R)
  y_pre <- y_obs - diag(Rinv)^(-1) * Rinv %*% y_obs
  return(y_pre)
}

getFitness <- function(model, sIn.vl = NULL, fIn.vl = NULL, sOut.vl) {
  # identify required statistic based on data
  stat <- "Q2cv" # Q2ext
  switch(stat,
         "Q2cv" = {# 1: leave-one-out cross-validation Q2
           y.hat <- getOut_loocv(model)
           eta <- 1 - (mean((model@sOut - y.hat)^2)/mean((model@sOut - mean(model@sOut))^2))
         },

         "Q2ext" = {# 3: External validation set Q2
           y.hat <- gaussian_cor(Ms, thetas) # correct here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           eta <- 1 - (mean((sOut.vl - y.hat)^2)/mean((sOut.vl - mean(model@sOut))^2))
         })

  return(eta)
}

getArgs <- function(ant, ds, df, sIn, fIn, base) {
  # recover base components
  s.state <- base$s.state
  f.state <- base$f.state
  f.dims <- base$f.dims
  f.dist <- base$f.dist
  f.fam <- base$f.fam
  k.type <- base$k.type

  # remove inactive scalar variables
  s.active.ac <- which(ant[1:ds] == 1) # index of active scalar inputs
  if (length(s.active.ac) > 0) {
    # print(paste("s.active:", s.active.ac)) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Remove this line
    sIn.ac <- sIn[,s.active.ac, drop = F]
  } else {
    sIn.ac <- NULL
  }
  piv <- ds

  defargs <- formals(funGp) # in case they are required below

  if (df > 0) {
    f.active.ac <- c() # index of active functional inputs
    f.dist.ac <- rep(0, df) # distance for each active functional input
    f.dims.ac <- rep(0, df) # dimension of each active functional input
    f.fam.ac <- rep(0, df) # basis function for each active functional input
    for (j in 1:df) {
      piv <- piv + 1
      if (ant[piv] == 1) f.active.ac <- c(f.active.ac, j)
      piv <- piv + 1
      if (ant[piv] == -1) { # assign default value
        f.dist.ac <- defargs$f_disType
        f.dims.ac <- defargs$f_pdims
        f.fam.ac <- defargs$f_family
        piv <- piv + 2
      } else {
        f.dist.ac[[j]] <- f.dist[[j]][ant[piv]]
        piv <- piv + 1
        # if (ant[piv] == 1) f.dims.ac[[j]] <- 0 else f.dims.ac[[j]] <- ant[piv] - 1
        f.dims.ac[[j]] <- ant[piv]
        piv <- piv + 1
        f.fam.ac[[j]] <- f.fam[[j]][ant[piv]]
      }
    }
    if (length(f.active.ac) > 0) {
      fIn.ac <- fIn[f.active.ac]
      f.dist.ac <- f.dist.ac[f.active.ac]
      f.dims.ac <- f.dims.ac[f.active.ac]
      f.fam.ac <- f.fam.ac[f.active.ac]

    } else { # set by default values
      fIn.ac <- defargs$fIn
      f.dist.ac <- defargs$f_disType
      f.dims.ac <- defargs$f_pdims
      f.fam.ac <- defargs$f_family
    }
  } else { # set by default values
    fIn.ac <- defargs$fIn
    f.dist.ac <- defargs$f_disType
    f.dims.ac <- defargs$f_pdims
    f.fam.ac <- defargs$f_family
    piv <- piv + df * 4
  }
  piv <- piv + 1
  k.type.ac <- k.type[ant[piv]]

  return(list(sIn = sIn.ac, fIn = fIn.ac, kerType = k.type.ac, f_disType = f.dist.ac, f_pdims = f.dims.ac, f_family = f.fam.ac))
}

# getArgs <- function() {
#   # remove inactive scalar variables
#   s.active.ac <- which(ants[i,1:ds] == 1) # index of active scalar inputs
#   if (length(s.active.ac) > 0) {
#     print(paste("s.active:", s.active.ac))
#     sIn.ac <- sIn[,s.active.ac, drop = F]
#   } else {
#     sIn.ac <- NULL
#   }
#   piv <- ds
#   f.active.ac <- c() # index of active functional inputs
#   f.dist.ac <- rep(0, df) # distance for each functional input
#   f.dims.ac <- rep(0, df) # dimension of each functional input
#   f.fam.ac <- rep(0, df) # basis function for each functional input
#   for (j in 1:df) {
#     piv <- piv + 1
#     if (ants[i, piv] == 1) f.active.ac <- c(f.active.ac, j)
#     piv <- piv + 1
#     f.dist.ac[[j]] <- f.dist[[j]][ants[i, piv]]
#     piv <- piv + 1
#     if (ants[i, piv] == 1) f.dims.ac[[j]] <- 0 else f.dims.ac[[j]] <- ants[i, piv] - 1
#     piv <- piv + 1
#     f.fam.ac[[j]] <- f.fam[[j]][ants[i, piv]]
#   }
#   if (length(f.active.ac) > 0) {
#     print(paste("f.active:", f.active.ac))
#     fIn.ac <- fIn[f.active.ac]
#     f.dist.ac <- f.dist.ac[f.active.ac]
#     f.dims.ac <- f.dims.ac[f.active.ac]
#     f.fam.ac <- f.fam.ac[f.active.ac]
#   } else {
#     fIn.ac <- NULL
#   }
#   piv <- piv + 1
#   k.type.ac <- k.type[ants[i, piv]]
#
#   return(list(sIn = sIn.ac, fIn = fIn.ac, kerType = k.type.ac, f_disType = f.dist.ac, f_pdims = f.dims.ac, f_family = f.fam.ac))
# }



sample.vec <- function(x, ...){
  x[sample(length(x), ...)]
}



updateSpace <- function(base = base, copy = copy, ds, df,
                        s_fixed = s_fixed, f_fixed = f_fixed, f_fixdim = f_fixdim, f_maxdim = f_maxdim,
                        f_setdist = f_setdist, f_setfam = f_setfam, k_settype = k_settype) {
  # recover originals
  s.state.0 <- base$s.state
  f.state.0 <- base$f.state
  f.dims.0 <- base$f.dims
  f.dist.0 <- base$f.dist
  f.fam.0 <- base$f.fam
  k.type.0 <- base$k.type

  # recover copies
  s.state.1 <- copy$s.state
  f.state.1 <- copy$f.state
  f.dims.1 <- copy$f.dims
  f.dist.1 <- copy$f.dist
  f.fam.1 <- copy$f.fam
  k.type.1 <- copy$k.type

  # update state of scalar inputs
  if (!is.null(s_fixed)) s.state.1[s_fixed] <- 1

  # update state of functional inputs
  if (!is.null(f_fixed)) f.state.1[f_fixed] <- 1

  # update the set of potential dimensions for functional inputs
  if (!is.null(f_fixdim)) {
    for (i in ncol(f_fixdim)) {
      f.dims.1[[f_fixdim[1,i]]] <- f_fixdim[2,i]
    }
  }

  # update the set of maximum dimensions
  if (!is.null(f_maxdim)) {
    for (i in ncol(f_maxdim)) {
      f.dims.1[[f_maxdim[1,i]]] <- 1:f_maxdim[2,i]
    }
  }

  # update the set of potential distances for functional inputs
  if (!is.null(f_setdist)) {
    ids <- as.numeric(names(f_setdist))
    for (i in length(ids)) {
      f.dist.1[[ids[i]]] <- f_setdist[[i]]
    }
  }

  # update the set of potential bases for functional inputs
  if (!is.null(f_setfam)) {
    ids <- as.numeric(names(f_setfam))
    for (i in length(ids)) {
      f.fam.1[[ids[i]]] <- f_setfam[[i]]
    }
  }

  # update the set of potential kernel functions
  if (!is.null(k_settype)) {
    k.type.1 <- k_settype
  }

  # return updated space
  clean <- list(s.state.1 = s.state.1, f.state.1 = f.state.1, f.dims.1 = f.dims.1,
                f.dist.1 = f.dist.1, f.fam.1 = f.fam.1, k.type.1 = k.type.1)

  return(clean)
}


checkAndPrint <- function(ds, df, space) {
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
