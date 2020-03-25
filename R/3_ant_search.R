# Method to plot a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @importFrom scales alpha
#' @importFrom graphics axis
#' @importFrom stats median setNames
#' @importFrom utils txtProgressBar setTxtProgressBar
# run_ACO <- function(sIn, fIn, sOut, ind.vl, param, env, base, extargs, time.str, time.lim, quietly, par.clust) {
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


  # recover pheromones and visibility
  # phero <- env$phero
  # visib <- env$visib

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
  # visib <- lapply(visib, function(M) M*0) # ojooooooooooooooooooooooooooo esto es temporal!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # browser()
  # run the colony
  for (c.gen in 1:n.iter) {
    # start_time <- Sys.time()
    cat(paste("Dispatching colony", c.gen, "\n"))
    # pb <- txtProgressBar(min = 0, max = n.pop, style = 3)

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

      # print(paste("id.ant:", id.ant))
      # print(myant$sol)

      # add step to selected ant based on decision rules
      # print("----------------< Next step")
      if (runif(1) <= q0) {
        # print("greedy")
        # antup <- nextNode_ACO(myant, "greedy", phero, visib, alp, bet, c.gen)
        antup <- nextNode_ACO(myant, "greedy", phero, c.gen)
      } else {
        # print("propor")
        # antup <- nextNode_ACO(myant, "propor", phero, visib, alp, bet, c.gen)
        antup <- nextNode_ACO(myant, "propor", phero, c.gen)
      }

      # print(paste("id.ant:", id.ant))
      # print(antup$sol)

      # perform local pheromone update only if the factor was not already fixed
      # print("----------------< Local update")
      # if (antup$layer - myant$layer == 1) {
      phero <- localUpd_ACO(phero, myant, antup, rho.l, tao0) #dt.l)
      # }

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
# browser() ##### estamos haciend pruebas aca!!
    # compute fitness of each ant
    fitness <- rep(0, n.pop)

    # browser()
    if (is.null(ind.vl)) {
      # cl <- parallel::makeCluster(3)
      res <- eval_loocv_ACO(sIn, fIn, sOut, extargs, base, ants, time.str, time.lim, quietly, par.clust)
      # parallel::stopCluster(cl)
    } else {
      # par.clust <- parallel::makeCluster(3)
      res <- eval_houtv_ACO(sIn, fIn, sOut, extargs, base, ants, ind.vl, time.str, time.lim, quietly, par.clust)
      # parallel::stopCluster(par.clust)
    }
    # browser()

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

    # browser()
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
    # browser()
    evol.fitness[c.gen] <- b.fitness
    all.ants[[c.gen]] <- ants
    all.fitness[[c.gen]] <- fitness

    dt <- difftime(Sys.time(), time.str, units = 'secs')
    if (dt >= time.lim) {
      cat(paste("\n** Time limit reached, exploration stopped after", format(as.numeric(dt), digits = 3, nsmall = 2), "seconds.\n"))
      break
    }
  }
# browser()

  # buena grafica --------------------<
  # plot(1, type = "n", xlab = "Colony", ylab = "Fitness", xlim = c(1, (c.gen + .3)), ylim = c(0, 1), xaxt = "n")
  # # axis(1, 1:c.gen)
  # axis(1, axtags(c.gen))
  # a <- rep(0, c.gen)
  # m <- 0
  # for (i in 1:c.gen) {
  #   points(rep(i, length(all.fitness[[i]])), all.fitness[[i]], pch = 21, bg = alpha("red", .4), col = alpha("red", .4))
  #   points(i, median(all.fitness[[i]]), pch = 21, bg = "blue", col = alpha("blue", .4))
  #   a[i] <- max(c(m, all.fitness[[i]]))
  #   m <- a[i]
  #   points(i, a[i], pch = 21, cex = 2, bg = NA, col = "magenta")
  # }
  # lines(a, col = "magenta")
  # buena grafica --------------------<

  # vrs <- sapply(all.ants, function(M) sum(apply(M, 2, var)))
  # plot(vrs, type = "b", ylim = c(0, max(vrs)))
  # abline(h = 0, col = "magenta")
  # browser()
  # buena grafica --------------------<
  # totDif <- sapply(all.ants, function(M) sum(apply(M, 2, function(v) length(unique(v))-1)))
  # plot(totDif, type = "b", ylim = c(0, max(totDif)))
  # abline(h = 0, col = "magenta")
  # buena grafica --------------------<
  # browser()

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

# nextNode_ACO <- function(myant, rule, phero, visib, alp, bet, c.gen) {
nextNode_ACO <- function(myant, rule, phero, c.gen) {
  # browser()
  # copy my ant to make updates on it
  antup <- myant

  # compute attractiveness
  layer <- myant$layer + 1
  level <- myant$level
  # tao <- unname(phero[[layer]][level,]) # level of pheromones of neighbors
  # nab <- unname(visib[[layer]][level,]) # visibility  of neighbors
  # attr <- alp * tao + bet * nab # attractiveness of each neighbor node
  attr <- unname(phero[[layer]][level,]) # level of pheromones of neighbors

  # if (grepl("Dist", names(phero)[layer])) { # ojoooooooooooooooooooooooooo eliminar!!!!!!!!!!!!!!!!!!!!!!!!!
  #   browser()
  # }
  # if (grepl("Dim", names(phero)[layer])) { # ojoooooooooooooooooooooooooo eliminar!!!!!!!!!!!!!!!!!!!!!!!!!
  #   attr <- attr[1:6]
  #   attr[1] <- 0
  # }

  # if (any(attr < 0)) {
  #   browser()
  # }

  # if (grepl("Dim", names(phero)[layer])) { # ojoooooooooooooooooooooooooo eliminar!!!!!!!!!!!!!!!!!!!!!!!!!
  #   browser()
  # }

  # if (grepl("State X1", names(phero)[layer])) { # ojoooooooooooooooooooooooooo eliminar!!!!!!!!!!!!!!!!!!!!!!!!!
  #   browser()
  # }

  # select next step
  switch(rule,
         "greedy" = {# choose the neighbor with largest attractiveness
           sel.level <- sample.vec(which(attr == max(attr)), 1)
         },

         "propor" = {# choose the neighbor randomly based on attractiveness pie
           sel.level <- sample(1:length(attr), 1, prob = attr/sum(attr))
         })

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
    antup$sol[x.skip] <- -1 # fill solution at skipped layers (arbitrarily 1 but has no effect)
    antup$layer <- tar.lay # update layer according to skip
    antup$level <- 1 # locate the ant arbitrarily at level 1 (irrelevant for next decision)
  }

  # return the updated ant
  return(antup)
}

localUpd_ACO <- function(phero, myant, antup, rho.l, tao0) { #dt.l) {
  # browser()
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
          # print(paste("phero0:", phero[[antup$layer]][o,d])) # flaggggggggggggggggggggggggggggggggggggg
          phero[[antup$layer]][o,d] <- (1 - rho.l) * phero[[antup$layer]][o,d] + rho.l * tao0 #(dt.l)
          # print(paste("phero1:", phero[[antup$layer]][o,d])) # flaggggggggggggggggggggggggggggggggggggg
        }
      } else {
        # print(paste("phero0:", phero[[antup$layer]][o,d])) # flaggggggggggggggggggggggggggggggggggggg
        phero[[antup$layer]][o,d] <- (1 - rho.l) * phero[[antup$layer]][o,d] + rho.l * tao0 #(dt.l)
        # print(paste("phero1:", phero[[antup$layer]][o,d])) # flaggggggggggggggggggggggggggggggggggggg
      }
    } else {
      # print(paste("phero0:", phero[[antup$layer]][o,d])) # flaggggggggggggggggggggggggggggggggggggg
      phero[[antup$layer]][o,d] <- (1 - rho.l) * phero[[antup$layer]][o,d] + rho.l * tao0 #(dt.l)
      # print(paste("phero1:", phero[[antup$layer]][o,d])) # flaggggggggggggggggggggggggggggggggggggg
    }
  } else {
    # print(paste("fixed:", phero[[antup$layer]][o,d])) # flaggggggggggggggggggggggggggggggggggggg
  }
  return(phero)
}

getElite_ACO <- function(fitness, n.ibest, ants, u.gbest, c.gen, b.ant, b.fitness){
  # browser()
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

globalUpd_ACO <- function(b.ants, b.fitness, phero, rho.g, tao0) {
  # browser()
  n.best <- nrow(b.ants)
  for (i in 1:n.best) {
    c.ant <- b.ants[i,] # current ant (dynamic during the outer loop)
    c.lev <- 1 # current level (dynamic during the inner loop)
    c.lay <- 1 # current layer (dynamic during the inner loop)
    # for (c.lay in 1:length(phero)) { # current layer (dynamic during the inner loop)
    while (c.lay <= length(phero)) {
      o <- c.lev
      if (grepl("Dim", names(phero)[c.lay])) {
        d <- c.ant[c.lay] + 1 # to correct projection dimension
      } else {
        d <- c.ant[c.lay]
      }

      # if (c.lev < 1) {
      #   browser()
      # }
      if (c.ant[c.lay] > 0) {
        if (phero[[c.lay]][o,d] < 1) {
          if (all(grepl("State F", names(phero)[c.lay]), d == 2)) {
            phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
            # identify the next layer of interest (just before next selection)
            # c.lay <- (which(!grepl("FI", names(phero))) > c.lay)[1]
            # browser()
            nff.lay <- which(!grepl("FI", names(phero)))
            # identify the next layer of interest (just before next selection)
            c.lay <- (nff.lay[nff.lay > c.lay])[1] - 1
            c.lev <- 1 # current level update

          } else {
            if (grepl("State X", names(phero)[c.lay])) { # if a scalar input has been activated
              # if (d == 1) {
              phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
              # }
            } else if (grepl("State F", names(phero)[c.lay])) { # if a functional input has been activated
              # if (d == 1) {
              if (c.ant[(c.lay-1)] == -1) { # if we come from an inactive functional input
                for (o in 1:nrow(phero[[c.lay]])) {
                  phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
                }
              } else {
                phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
              }
              # }
            } else { # current level does not have anything to do with the state of an input
              if (c.ant[(c.lay-1)] == -1) { # if we come from an inactive functional input
                for (o in 1:nrow(phero[[c.lay]])) {
                  phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
                }
              } else {
                phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],tao0) # pheromone update
              }
            }
            # if (c.ant[(c.lay-1)] == -1) {
            #   for (o in 1:nrow(phero[c.lay-1])) {
            #     phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],1e-8) # pheromone update
            #   }
            # } else if (all(!grepl("State X", names(phero)[c.lay]), d == 1)) {
            #   phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],1e-8) # pheromone update
            # }
            # print(paste("phero0:", phero[[c.lay]][o,d]))
            # phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * max(b.fitness[i],1e-8) # pheromone update
            c.lev <- d # current level update
            # print(paste("phero1:", phero[[c.lay]][o,d]))
          }
        } else {
          # print(paste("fixed:", phero[[c.lay]][o,d]))
          # browser()
          c.lev <- d
        }
      }
      c.lay <- c.lay + 1
      # c.tao <- phero[[c.lay]][o,d] # current pheromone value
      # p.tao <- (1 - rho.g) * c.tao + rho.g * b.fitness[i] # potential tao # p.tao <- (1 - rho.g) * c.tao + (1 - (1 - rho.g) * c.tao) * b.fitness[i] # potential tao
      # phero[[c.lay]][o,d] <- max(c.tao, p.tao) # pheromone update
    }
  }
  return(phero)
}
