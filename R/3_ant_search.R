# ==========================================================================================================
# Skeleton of the Ant Colony Optimization algoritm for model selection
# ==========================================================================================================
run_ACO <- function(sIn, fIn, sOut, ind.vl, param, phero, base, extargs, time.str, time.lim, quietly, par.clust) {
  # recover heuristic parameters
  #___________________________________________________________________________________________
  # <---> population factors
  n.iter <- param$n.iter
  n.pop <- param$n.pop

  # <---> transition rules
  q0 <- param$q0
  # alp <- param$alp
  # bet <- param$bet

  # <---> local pheromone update
  rho.l <- param$rho.l
  # dt.l <- param$dt.l

  # <---> global pheromone update
  tao0 <- param$tao0
  u.gbest <- param$u.gbest
  n.ibest <- param$n.ibest
  rho.g <- param$rho.g
  #___________________________________________________________________________________________

  # set up heuristic controllers and statistics
  #___________________________________________________________________________________________
  # <---> required controllers
  n.layers <- length(phero) # number of layers in the network/factors in the experiment (some may be fixed constant)
  timestop <- F

  # <---> for statistics
  evol.fitness <- rep(0, n.iter) # fitness of best ant of each colony
  all.ants <- list() # log of all explored ants
  # all.fitness <- matrix(nrow = n.pop, ncol = n.iter) # fitness of all explored ants
  all.fitness <- list() # fitness of all explored ants
  crashes <- list()

  # <---> best solution and its fitness
  b.ant <- "just use the mean" # initialize best ant
  b.fitness <- 0 # initialize best fitness
  #___________________________________________________________________________________________

  # run the colony
  for (c.gen in 1:n.iter) {
    # start_time <- Sys.time()
    cat(paste("Dispatching colony", c.gen, "\n"))

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

      # add step to selected ant based on decision rules
      if (runif(1) <= q0) {
        antup <- nextNode_ACO(myant, "greedy", phero, c.gen)
      } else {
        antup <- nextNode_ACO(myant, "propor", phero, c.gen)
      }

      # perform local pheromone update only if the factor was not already fixed
      phero <- localUpd_ACO(phero, myant, antup, rho.l, tao0)

      # save the ant
      ants[antup$id,] <- antup$sol
      antsLayers[antup$id] <- antup$layer
      antsLevels[antup$id] <- antup$level

      # check if the ant is done
      if (antup$layer == n.layers) {
        active <- getActiveIn_ACO(antup$sol, sIn, fIn, base)
        if ((length(active$s.active) + length(active$f.active)) > 0) {
          # tag ant as complete
          antsPartial[antup$id] <- F

        } else {
          # restart ant
          ants[antup$id,] <- NA
          antsLayers[antup$id] <- 0
          antsLevels[antup$id] <- 1
        }
      }

      # check if the colony is done
      if (all(!antsPartial)) colPartial <- F
    }

    # compute fitness of each ant
    fitness <- rep(0, n.pop)
    if (is.null(ind.vl)) {
      res <- eval_loocv_ACO(sIn, fIn, sOut, extargs, base, ants, time.str, time.lim, quietly, par.clust)
    } else {
      res <- eval_houtv_ACO(sIn, fIn, sOut, extargs, base, ants, ind.vl, time.str, time.lim, quietly, par.clust)
    }

    # extract complete evaluations
    done <- which(!sapply(res, is.null))
    argsList <- lapply(res[done], `[[`, 1)
    modelList <- lapply(res[done], `[[`, 2)
    fitness <- sapply(res[done], `[[`, 3)

    # identify crashes and usable models
    ids.cr <- which(is.na(fitness))
    ids.ok <- which(!is.na(fitness))

    # save args of crashes in ant mode (if any)
    if (length(ids.cr) > 0) {
      crashes[[c.gen]] <- ants[ids.cr,,drop = F]
    }

    # extract ants and fitness for global update
    fitness <- fitness[ids.ok]
    ants <- ants[ids.ok,]
    elite <- getElite_ACO(fitness, n.ibest, ants, u.gbest, c.gen, b.ant, b.fitness)

    # perform global pheromone update
    phero <- globalUpd_ACO(elite$ants.up, elite$fitness.up, phero, rho.g, tao0)

    # save best ant, fitness, agrs and model
    if (fitness[elite$b.ind[1]] > b.fitness) {
      b.ant <- ants[elite$b.ind[1],]
      b.fitness <- fitness[elite$b.ind[1]]
      b.args <- argsList[ids.ok[elite$b.ind[1]]][[1]]
      b.model <- modelList[ids.ok[elite$b.ind[1]]][[1]]
    }

    # save data for statistics
    evol.fitness[c.gen] <- b.fitness
    all.ants[[c.gen]] <- ants
    all.fitness[[c.gen]] <- fitness

    dt <- difftime(Sys.time(), time.str, units = 'secs')
    if (dt >= time.lim) {
      cat(paste("\n** Time limit reached, exploration stopped after", format(as.numeric(dt), digits = 3, nsmall = 2), "seconds.\n"))
      break
    }
  }

  # merge all successful ants
  all.ants <- do.call(rbind, all.ants)

  # identify duplicates (if any) and keep only the best
  all.fitness.v <- unlist(all.fitness)
  ord.fitness <- sort(all.fitness.v, decreasing = T)
  ord.ants <- all.ants[order(all.fitness.v, decreasing = T),]
  top.fitness <- ord.fitness[!duplicated(ord.ants)]
  top.ants <- ord.ants[!duplicated(ord.ants),]
  top.fitness <- sort(top.fitness, decreasing = T)
  top.ants <- top.ants[order(top.fitness, decreasing = T),]

  # remove duplicates in crashed ants
  if (length(crashes) > 0) {
    crashes <- do.call(rbind, crashes)
    crashes <- crashes[!duplicated(crashes),,drop = F]
  }

  cat("\nAnts are done ;)\n")

  return(list(model = b.model, sol.vec = b.ant, sol.args = b.args, b.fitness = b.fitness,
              log.suc = top.ants, log.fitness = top.fitness, log.cra = crashes,
              all.details = list(param = param, evolution = all.fitness)))
}
# ==========================================================================================================



# ==========================================================================================================
# Function for the selection of next node in the decision network for ants
# ==========================================================================================================
nextNode_ACO <- function(myant, rule, phero, c.gen) {
  # copy my ant to make updates on it
  antup <- myant

  # compute attractiveness
  layer <- myant$layer + 1
  level <- myant$level
  attr <- unname(phero[[layer]][level,]) # level of pheromones of neighbors

  # select next step
  switch(rule,
         "greedy" = {# choose the neighbor with largest attractiveness
           sel.level <- sample.vec(which(attr == max(attr)), 1)
         },

         "propor" = {# choose the neighbor randomly based on attractiveness pie
           sel.level <- sample(1:length(attr), 1, prob = attr/sum(attr))
         })

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
    nff.lay <- which(!grepl("FI", names(phero)))
    tar.lay <- (nff.lay[nff.lay > layer])[1] - 1
    x.skip <- (layer+1):tar.lay # get index of the layers to skip
    antup$sol[x.skip] <- -1 # fill solution at skipped layers (arbitrarily 1 but has no effect)
    antup$layer <- tar.lay # update layer according to skip
    antup$level <- 1 # locate the ant arbitrarily at level 1 (irrelevant for next decision)
  }

  # return the updated ant
  return(antup)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to perform local pheromone update
# ==========================================================================================================
localUpd_ACO <- function(phero, myant, antup, rho.l, tao0) {
  o <- myant$level
  if (antup$layer - myant$layer == 1) {
    d <- antup$level
  } else {
    d <- 2
    antup$layer <- myant$layer + 1
  }

  if (phero[[antup$layer]][o,d] < 1) {
    if (antup$layer > 1) {
      if (antup$sol[antup$layer-1] == -1) {
        for (o in 1:nrow(phero[[antup$layer]])) {
          phero[[antup$layer]][o,d] <- (1 - rho.l) * phero[[antup$layer]][o,d] + rho.l * tao0
        }
      } else {
        phero[[antup$layer]][o,d] <- (1 - rho.l) * phero[[antup$layer]][o,d] + rho.l * tao0
      }
    } else {
      phero[[antup$layer]][o,d] <- (1 - rho.l) * phero[[antup$layer]][o,d] + rho.l * tao0
    }
  } else {
  }
  return(phero)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to identify best ants for global pheromone update
# ==========================================================================================================
getElite_ACO <- function(fitness, n.ibest, ants, u.gbest, c.gen, b.ant, b.fitness){
  # identify best n.ibest ants
  b.ind <- order(fitness, decreasing = T)[1:min(n.ibest, length(fitness))]

  # remove duplicates if there is any
  if (n.ibest > 1) {
    u.ind <- b.ind[!duplicated(ants[b.ind,])] # unique best ants
  } else {
    u.ind <- b.ind
  }

  # group best ants and their fitness
  if (!is.matrix(ants)) ants <- matrix(ants, ncol = length(ants))
  ants.up <- ants[u.ind,,drop = F]
  fitness.up <- fitness[u.ind]

  # add the global best only if required
  if (all(u.gbest, c.gen > 1)) {
    ants.up <- rbind(ants.up, b.ant)
    fitness.up <- c(fitness.up, b.fitness)
  }

  return(list(ants.up = ants.up, fitness.up = fitness.up, b.ind = b.ind))
}
# ==========================================================================================================



# ==========================================================================================================
# Function to perform global pheromone update
# ==========================================================================================================
globalUpd_ACO <- function(b.ants, b.fitness, phero, rho.g, tao0) {
  n.best <- nrow(b.ants)
  for (i in 1:n.best) {
    c.ant <- b.ants[i,] # current ant (dynamic during the outer loop)
    c.lev <- 1 # current level (dynamic during the inner loop)
    c.lay <- 1 # current layer (dynamic during the inner loop)
    while (c.lay <= length(phero)) {
      o <- c.lev
      if (grepl("Dim", names(phero)[c.lay])) {
        d <- c.ant[c.lay] + 1 # to correct projection dimension
      } else {
        d <- c.ant[c.lay]
      }

      if (c.ant[c.lay] > 0) {
        if (phero[[c.lay]][o,d] < 1) {
          if (all(grepl("State F", names(phero)[c.lay]), d == 2)) {
            phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
            # identify the next layer of interest (just before next selection)
            nff.lay <- which(!grepl("FI", names(phero)))
            c.lay <- (nff.lay[nff.lay > c.lay])[1] - 1
            c.lev <- 1 # current level update

          } else {
            if (grepl("State X", names(phero)[c.lay])) { # if a scalar input has been activated
              phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update

            } else if (grepl("State F", names(phero)[c.lay])) { # if a functional input has been activated
              if (c.ant[(c.lay-1)] == -1) { # if we come from an inactive functional input
                for (o in 1:nrow(phero[[c.lay]])) {
                  phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
                }

              } else {
                phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
              }

            } else { # current level does not have anything to do with the state of an input
              if (c.ant[(c.lay-1)] == -1) { # if we come from an inactive functional input
                for (o in 1:nrow(phero[[c.lay]])) {
                  phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
                }

              } else {
                phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
              }
            }
            c.lev <- d # current level update
          }
        } else {
          c.lev <- d
        }
      }
      c.lay <- c.lay + 1
    }
  }
  return(phero)
}
# ==========================================================================================================
