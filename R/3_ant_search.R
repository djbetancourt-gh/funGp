# ==========================================================================================================
# Ants log
# ==========================================================================================================
#' @title Class: data structures related to the kernel of a funGp model
#' @description Fill this!!!!!!!!!
#'
#' @slot sols Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#' @slot args Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#' @slot fitness Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setClass("antsLog",
         representation(
           sols = "data.frame",            # kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}
           args = "list",                  # distance type. To be chosen from {"scalar", "functional"}
           fitness = "numeric"             # search method
         ),
         validity = function(object) {T})
# ==========================================================================================================


# Method to plot a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @importFrom scales alpha
#' @importFrom graphics axis
#' @importFrom stats median setNames
#' @importFrom utils setTxtProgressBar txtProgressBar
run_ACO <- function(sIn, fIn, sOut, ind.vl, param, env, base, extargs) {
  # recover heuristic parameters
  #___________________________________________________________________________________________
  # <---> population factors
  n.gen <- param$n.gen
  n.pop <- param$n.pop

  # <---> transition rules
  q0 <- param$q0
  alp <- param$alp
  bet <- param$bet

  # <---> local pheromone update
  rho.l <- param$rho.l
  dt.l <- param$dt.l

  # <---> global pheromone update
  u.gbest <- param$u.gbest
  n.lbest <- param$n.lbest
  rho.g <- param$rho.g
  #___________________________________________________________________________________________


  # recover pheromones and visibility
  phero <- env$phero
  visib <- env$visib


  # set up heuristic controllers and statistics
  #___________________________________________________________________________________________
  # <---> required controllers
  n.layers <- length(phero) # number of layers in the network/factors in the experiment (some may be fixed constant)

  # <---> for statistics
  evol.fitness <- rep(0, n.gen) # fitness of best ant of each colony
  all.ants <- list() # log of all explored ants
  all.fitness <- matrix(nrow = n.pop, ncol = n.gen) # fitness of all explored ants

  # <---> best solution and its fitness
  b.ant <- "just use the mean" # initialize best ant
  b.fitness <- 0 # initialize best fitness
  #___________________________________________________________________________________________

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
        antup <- nextNode_ACO(myant, "greedy", phero, visib, alp, bet, c.gen)
      } else {
        # print("propor")
        antup <- nextNode_ACO(myant, "propor", phero, visib, alp, bet, c.gen)
      }

      # if (antup$layer == 8) {
      #   browser()
      #   if (antup$sol[8] != 2) {
      #     browser()
      #   }
      # }


      # perform local pheromone update only if the factor was not already fixed
      # print("----------------< Local update")
      if (antup$layer - myant$layer == 1) {
        phero <- localUpd_ACO(phero, myant, antup, rho.l, dt.l)
      }

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

    # if there are no active inputs left, restart the ant !!!!!!!!!!!!!!!!!!!!!!! complete this!!!!!!!!!!!!!!!!!!!!!!!!!
    # if (all(is.null(args$sIn), is.null(args$fIn))) browser()
    # getActiveIn_ACO(ants[i,], sIn, fIn, base)

    # compute fitness of each ant
    fitness <- rep(0, n.pop)
    for (i in 1:n.pop) {
      if (is.null(ind.vl)) {
        # translate ant data into funGp arguments format
        args <- formatSol_ACO(ants[i,], sIn, fIn, base)

        # build the model
        model <- quiet(funGp(sIn = args$sIn, fIn = args$fIn, sOut = sOut, kerType = args$kerType,
                             f_disType = args$f_disType, f_pdims = args$f_pdims, f_basType = args$f_basType,
                             nugget = extargs$nugget, n.starts = extargs$n.starts, n.presample = extargs$n.presample))

        # compute model fitness
        fitness[i] <- max(getFitness(model),0)
      } else {
        n.rep <- ncol(ind.vl)# number of replicates
        rep.fitness <- rep(0, n.rep)
        for (r in 1:n.rep) {
          # split data into training and validation
          data <- splitData(sIn, fIn, sOut, ind.vl[,r]) # aca no hay que poner entradas inactivas en null

          # translate ant data into funGp arguments format
          args <- formatSol_ACO(ants[i,], sIn = data$sIn.tr, fIn = data$fIn.tr, base) # esto se encarga de poner en null las entradas inact

          # build the model
          model <- quiet(funGp(sIn = args$sIn, fIn = args$fIn, sOut = data$sOut.tr, kerType = args$kerType,
                               f_disType = args$f_disType, f_pdims = args$f_pdims, f_basType = args$f_basType,
                               nugget = extargs$nugget, n.starts = extargs$n.starts, n.presample = extargs$n.presample))

          # identify active inputs of both types
          active <- getActiveIn_ACO(ants[i,], sIn, fIn, base)

          # compute model-replicate fitness
          rep.fitness[r] <- max(getFitness(model, data$sIn.vl, data$fIn.vl, data$sOut.vl, active),0)
        }

        # compute average model fitness
        fitness[i] <- mean(rep.fitness)
      }

      # # translate ant data into funGp arguments format
      # args <- formatSol_ACO(ants[i,], sIn, fIn, base)
      #
      # if (all(is.null(args$sIn), is.null(args$fIn))) browser()
      #
      # # build the model
      # model <- quiet(funGp(sIn = args$sIn, fIn = args$fIn, sOut = sOut, kerType = args$kerType,
      #                      f_disType = args$f_disType, f_pdims = args$f_pdims, f_basType = args$f_basType))
      #
      # # compute model fitness #!!!!!!!!!!!!!!!!!! extend to other metrics
      # # getFitness(model, sIn.vl, fIn.vl, sOut.vl)
      # if (is.null(ind.vl)) {
      #   fitness[i] <- max(getFitness(model),0)
      # } else {
      #   fitness[i] <- max(getFitness(model, sIn, fIn, ind.vl, getActiveIn_ACO(ants[i,], sIn, fIn, base)),0)
      # }
      # # print(ants[i,])
      # # print(fitness[i])
      #
      # save the model if it is the global best
      if (fitness[i] > b.fitness) {
        b.model <- model
        b.args <- args
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)

    # extract ants and fitness for global update
    res <- getElite_ACO(fitness, n.lbest, ants, u.gbest, c.gen, b.ant, b.fitness)

    # perform global pheromone update
    # print("----------------< Global update")
    phero <- globalUpd_ACO(res$ants.up, res$fitness.up, phero, rho.g)

    # save best ant
    if (fitness[res$b.ind[1]] > b.fitness) {
      b.ant <- ants[res$b.ind[1],]
      b.fitness <- fitness[res$b.ind[1]]
    }
    # print(b.ant)

    print(Sys.time() - start_time)
    cat("\n")

    # save data for statistics
    evol.fitness[c.gen] <- b.fitness
    all.ants[[c.gen]] <- ants
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
    points(rep(i, n.pop), all.fitness[,i], pch = 21, bg = alpha("red", .4), col = alpha("red", .4))
    points(i, median(all.fitness[,i]), pch = 21, bg = "blue", col = alpha("blue", .4))
    # legend(x = (i + .2), y = (evol.fitness[i] - .05), legend = format(evol.fitness[i], digits = 2, nsmall = 3), cex = 1,
    #        xjust = 0.5,      # 0.5 means center adjusted
    #        yjust = 0.5,      # 0.5 means center adjusted
    #        x.intersp = -0.5, # adjust character interspacing as you like to effect box width
    #        y.intersp = 0.1,  # adjust character interspacing to effect box height
    #        adj = c(0, 0.5))
  }
  # lines(evol.fitness, lty = 2, col = "blue")
  points(evol.fitness, pch = 21, bg = alpha("blue", 0), col = "black")

  cat("\nAnts are done ;)")

  # merge all ants
  all.ants <- do.call(rbind, all.ants)

  # identify duplicates (if any) and keep only the best
  ord.fitness <- sort(all.fitness, decreasing = T)
  ord.ants <- all.ants[order(all.fitness, decreasing = T),]
  top.fitness <- ord.fitness[!duplicated(ord.ants)]
  top.ants <- ord.ants[!duplicated(ord.ants),]
  top.fitness <- sort(top.fitness, decreasing = T)
  top.ants <- top.ants[order(top.fitness, decreasing = T),]

  ####### not sure that b.ant matches b.args

  return(list(model = b.model, sol.vec = b.ant, sol.args = b.args, b.fitness = b.fitness,
              log.vec = top.ants, log.fitness = top.fitness, details = param))
}

nextNode_ACO <- function(myant, rule, phero, visib, alp, bet, c.gen) {
  # browser()
  # copy my ant to make updates on it
  antup <- myant

  # compute attractiveness
  layer <- myant$layer + 1
  level <- myant$level
  tao <- unname(phero[[layer]][level,]) # level of pheromones of neighbors
  nab <- unname(visib[[layer]][level,]) # visibility  of neighbors
  attr <- alp * tao + bet * nab # attractiveness of each neighbor node

  # select next step
  switch(rule,
         "greedy" = {# choose the neighbor with largest attractiveness
           sel.level <- sample.vec(which(attr == max(attr)), 1)
         },

         "propor" = {# choose the neighbor randomly based on attractiveness pie
           sel.level <- sample(1:length(tao), 1, prob = attr/sum(attr))
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

localUpd_ACO <- function(phero, myant, antup, rho.l, dt.l) {
  # browser()
  o <- myant$level
  d <- antup$level
  # if (any(compareNA(antup$sol, -1))) {
  #   browser()
  # }
  if (phero[[antup$layer]][o,d] < 1) {
    # print(paste("phero0:", phero[[antup$layer]][o,d]))
    phero[[antup$layer]][o,d] <- (1 - rho.l) * phero[[antup$layer]][o,d] + rho.l * (dt.l)
    # print(paste("phero1:", phero[[antup$layer]][o,d]))
  } else {
    # print(paste("fixed:", phero[[antup$layer]][o,d]))
  }
  return(phero)
}

getElite_ACO <- function(fitness, n.lbest, ants, u.gbest, c.gen, b.ant, b.fitness){
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

globalUpd_ACO <- function(b.ants, b.fitness, phero, rho.g) {
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

      # if (c.lev < 1) {
      #   browser()
      # }
      if (c.ant[c.lay] > 0) {
        if (phero[[c.lay]][o,d] < 1) {
          if (all(grepl("State F", names(phero)[c.lay]), c.lev == 2)) {
            # identify the next layer of interest (just before next selection)
            # c.lay <- (which(!grepl("FI", names(phero))) > c.lay)[1]
            # browser()
            nff.lay <- which(!grepl("FI", names(phero)))
            # identify the next layer of interest (just before next selection)
            c.lay <- (nff.lay[nff.lay > c.lay])[1]
            c.lev <- 1 # current level update
          } else {
            # print(paste("phero0:", phero[[c.lay]][o,d]))
            phero[[c.lay]][o,d] <- (1 - rho.g) * phero[[c.lay]][o,d] + rho.g * b.fitness[i] # pheromone update
            c.lev <- d # current level update
            # print(paste("phero1:", phero[[c.lay]][o,d]))
          }
        } else {
          # print(paste("fixed:", phero[[c.lay]][o,d]))
          # browser()
          c.lev <- d
        }
      }
      # c.tao <- phero[[c.lay]][o,d] # current pheromone value
      # p.tao <- (1 - rho.g) * c.tao + rho.g * b.fitness[i] # potential tao # p.tao <- (1 - rho.g) * c.tao + (1 - (1 - rho.g) * c.tao) * b.fitness[i] # potential tao
      # phero[[c.lay]][o,d] <- max(c.tao, p.tao) # pheromone update
    }
  }
  return(phero)
}
