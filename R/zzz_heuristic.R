# ==========================================================================================================
# Class for factories to produce structre-optimized funGp models
# ==========================================================================================================



# ==========================================================================================================
# Developer oriented methods
# ==========================================================================================================

# Constructor of the class
# ----------------------------------------------------------------------------------------------------------
# @title Class: Fill!!!!
# @description Fill this!!!!!!!!!
#
# @slot ppp Object of class \code{"character"}. Fill!!!!!!!!!!
#
# @rdname factory-class
#
# @author José Betancourt, François Bachoc and Thierry Klein
# @export
# setClass("funGpFactory",
#          representation(
#            ppp = "character"              # a ppp
#          ),
#          validity = function(object) {T})
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
symbolicfunGp2 <- function(sIn = NULL, fIn = NULL, sOut, tr.ind = NULL) {
  print("Hey, lets optimize that modal!")

  ########## factors
  # the list of variables to optimize is not asked to the user, it is extracted from the data structures provided
  # kerType <- c("gauss", "matern5_2")
  # f_family <- c("B-splines", "PCA")
  # f_disType <- c("L2_bygroup", "L2_byindex")
  # f_tryodim <- c(T, F)
  # f_maxpdim <- c(4, 8) # this should only be taken into account if the person asked for "L2_byindex" distance. should match number of f inputs
  # solspace <- list(kerType = kerType, f_family = f_family, f_disType = f_disType, f_tryodim = f_tryodim, f_maxpdim = f_maxpdim)

  ########## fixed levels
  # s_active <- NULL#c(2) # scalar input 2 should allways be in the model
  # f_active <- c(1) # functional input 1 should always be in the model
  # f_family <- matrix(c(1,"PCA"), nrow = 2)
  # keepfixed <- list(s_active = s_active, f_active = f_active, f_family = f_family)

  # All levels


  ##### default values
  s.always <- c(F, F) # as many elements as scalar inputs
  f.always <- c(F, F) # as many elements as functional inputs
  f.dim <- c(NA, NA) # as many elements as functional inputs
  f.dmax <- c(NA, NA) # as many elements as functional inputs

  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(fIn)) df <- length(fIn) else df <- 0
  # full experiment
  s.state.c <- s.state <- c(0,0) # 0 = free, 1 = fixed
  f.state.c <- f.state <- c(0,0) # 0 = free, 1 = fixed
  f1.dim <- 0:ncol(fIn[[1]])
  f2.dim <- 0:ncol(fIn[[2]])
  f.dim.c <- f.dim <- list(f1.dim, f2.dim)
  f1.dst <- c("L2_byindex", "L2_bygroup")
  f2.dst <- c("L2_byindex", "L2_bygroup")
  f.dst.c <- f.dst <- list(f1.dst, f2.dst)
  f1.fam <- c("B-splines", "PCA")
  f2.fam <- c("B-splines", "PCA")
  f.fam.c <- f.fam <- list(f1.fam, f2.fam)
  k.type.c <- k.type <- c("gaussian", "matern5_2", "matern3_2")

  # cat("Experiment summary______________")
  checkAndPrint(s.state, f.state, f.dim, f.dst, f.fam, k.type)
  # knitr::kable(cbind(fixe.fact, fixe.values), align = "c", col.names = c("Factor", "Levels"), caption = "Fixed factors")

  # user specifies changes
  s_fixed <- c(1) # keep X1 always active
  f_fxdim <- matrix(c(2,4), ncol = 1) # set dimension of F2 at 4
  f_mxdim <- matrix(c(1,5), ncol = 1) # set max dimension of F1 at 5
  f_dst <- list("2" = c("L2_byindex")) # only test L2_index distance for F2
  f_fam <- list("1" = c("B-splines")) # only B-splines projection for F1
  k_type <- c("matern5_2", "matern3_2") # only matern5_2 and matern3_2 kernels

  # interpreting
  # --------------------------------------
  # state of scalar inputs
  if (!is.null(s_fixed)) {
    s.state[s_fixed] <- 1
  }
  # state of functional inputs
  # if (!is.null(f_fixed)) {
  #   f.state[f_fixed] <- 1
  # }
  if (!is.null(f_fxdim)) {
    for (i in ncol(f_fxdim)) {
      f.dim[[f_fxdim[1,i]]] <- f_fxdim[2,i]
    }
  }
  if (!is.null(f_mxdim)) {
    for (i in ncol(f_mxdim)) {
      f.dim[[f_mxdim[1,i]]] <- 1:f_mxdim[2,i]
    }
  }
  if (!is.null(f_dst)) {
    ids <- as.numeric(names(f_dst))
    for (i in length(ids)) {
      f.dst[[ids[i]]] <- f_dst[[i]]
    }
  }
  if (!is.null(f_fam)) {
    ids <- as.numeric(names(f_fam))
    for (i in length(ids)) {
      f.fam[[ids[i]]] <- f_fam[[i]]
    }
  }
  if (!is.null(k_type)) {
    k.type <- k_type
  }

  checkAndPrint(s.state, f.state, f.dim, f.dst, f.fam, k.type)

  # create visibility and pheromones lists
  visib <- phero <- list()

  s.vis0 <- .7 # probability of activating any scalar input
  s.phr0 <- .5 # initial pheromones charge for the node where a scalar input is activated

  # set up initial visibilities and pheromones for scalar inputs
  for (i in 1:ds) {
    if (s.state[i] == 0) { # free, distribute load
      visib[[i]] <- matrix(rep(c(s.vis0, (1-s.vis0)), 2), nrow = 2, byrow = T)
      phero[[i]] <- matrix(rep(c(s.phr0, (1-s.phr0)), 2), nrow = 2, byrow = T)
    } else { # fixed, put all in first colum
      visib[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
      phero[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
    }
  }

  for (j in 1:df) {
    # state of functional input j
    i <- length(visib) + 1
    if (f.state[j] == 0) { # free, distribute load
      visib[[i]] <- matrix(rep(c(s.vis0, (1-s.vis0)), 2), nrow = 2, byrow = T)
      phero[[i]] <- matrix(rep(c(s.phr0, (1-s.phr0)), 2), nrow = 2, byrow = T)
    } else { # fixed, put all in first colum
      visib[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
      phero[[i]] <- matrix(rep(c(1, 0), 2), nrow = 2, byrow = T)
    }

    # distance for functional input j
    i <- length(visib) + 1
    nr <- nrow(visib[[i-1]])
    if (length(f.dst[[j]]) > 1) { # free, distribute load
      v <- rep(.5, length(f.dst.c[[j]]))
      visib[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
      phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
    } else { # fixed, put all in preferred colum
      v <- rep(0, length(f.dst.c[[j]]))
      v[which(f.dst.c[[j]] == f.dst[[j]])] <- 1
      visib[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
      phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
    }

    # dimension of functional input j
    b <- .5
    i <- length(visib) + 1
    nr <- nrow(visib[[i-1]])
    if (length(f.dim[[j]]) > 1) { # free, distribute load
      vind <- vgro <- rep(0, length(f.dim.c[[j]]))
      if (0 %in% f.dim[[j]]) {
        l <- exp(-b * c(max(f.dim[[j]]), f.dim[[j]][-1]))
      } else {
        l <- exp(-b * f.dim[[j]])
      }
      vind[(f.dim[[j]]+1)] <- l
      vgro[(f.dim[[j]]+1)] <- 1
      visib[[i]] <- rbind(vgro/sum(vgro), vind/sum(vind))
      v <- rep(.5, length(f.dim.c[[j]]))
      phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
    } else { # fixed, put all in preferred colum
      v <- rep(0, length(f.dim.c[[j]]))
      v[which(f.dim.c[[j]] == f.dim[[j]])-1] <- 1
      visib[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
      phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
    }

    # basis for functional input j
    i <- length(visib) + 1
    nr <- nrow(visib[[i-1]])
    if (length(f.fam[[j]]) > 1) { # free, distribute load
      v <- rep(.5, length(f.fam.c[[j]]))
      visib[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
      phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
    } else { # fixed, put all in preferred colum
      v <- rep(0, length(f.fam.c[[j]]))
      v[which(f.fam.c[[j]] == f.fam[[j]])] <- 1
      visib[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
      phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
    }
  }

  # kernel function
  i <- length(visib) + 1
  nr <- nrow(visib[[i-1]])
  if (length(k.type) == length(k.type.c)) { # free, distribute load
    v <- rep(.5, length(k.type.c))
    visib[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
    phero[[i]] <- matrix(rep(v, nr), nrow = nr, byrow = T)
  } else { # fixed, put all in preferred colum
    v <- rep(0, length(k.type.c))
    l <- unname(sapply(k.type, function(x) which(k.type.c == x)))
    v[l] <- 1
    visib[[i]] <- matrix(rep(v/sum(v), nr), nrow = nr, byrow = T)
    phero[[i]] <- matrix(rep(v/sum(v), nr), nrow = nr, byrow = T)
  }

  visib
  phero
  browser()
  # replace the rep by the number of original levels

  # hacer que cada fila siempre sume 1!!















  # user inputs
  s_always <- c(1)
  f_dim <- matrix(c(2,4), nrow = 2)
  f_dmax <- matrix(c(1,5), nrow = 2)

  # updates based on user inputs
  s.always[s_always] <- T
  f.dim[f_dim[1,]] <- f_dim[2,]
  f.dmax[f_dmax[1,]] <- f_dmax[2,]


  print("ok")










  # identity factors and levels, and create the required data structures
  # 1. active scalar inputs
  # 2. active functional inputs
  # 3. types of distances for functions
  # 4. basis family for functions
  # 5. projection dimension
  # 6. type of kernel function
  browser()

  # 1. set up everything about scalars
  # sIn: ds
  # solspace: --
  # keepfixed: active inputs
  if (!is.null(sIn)) ds <- ncol(sIn) else ds <- 0
  if (!is.null(keepfixed$s_active)) {
    s.indfix <- keepfixed$s_active
    s.indvar <- c(1:ds)[-s.indfix]
  } else {
    s.indvar <- c(1:ds)
  }
  s.nvar <- length(s.indvar)

  # 1. set up everything about functions
  # fIn: df, k
  # solspace: basis, distance, try original dim, max dim
  # keepfixed: active inputs, basis for functions
  if (!is.null(fIn)) df <- length(fIn) else df <- 0
  if (!is.null(keepfixed$f_active)) {
    f.indfix <- keepfixed$f_active
    f.indvar <- c(1:df)[-f.indfix]
  } else {
    f.indvar <- c(1:df)
  }
  f.nvar <- length(f.indvar)

  if (!is.null(keepfixed$f_pdim)) {

  }



  # create visibility and pheromones lists
  visib <- phero <- list()

  s.vis0 <- .7 # probability of activating any scalar input
  s.phr0 <- .5 # initial pheromones charge for the node where a scalar input is activated

  if (all(!is.null(sIn), !is.null(fIn))) { # Hybrid-input case *******************************************
    # extract information from user inputs specific to the hybrid-input case
    sIn <- as.matrix(sIn)
    ds <- ncol(sIn)
    df <- length(fIn)
    f_dims <- sapply(fIn, ncol)

    # set up initial visibilities and pheromones for scalar inputs
    if (1 %in% s.indvar) { # distribute load
      visib[[1]] <- matrix(c(s.vis0, (1-s.vis0)), nrow = 1)
      phero[[1]] <- matrix(c(s.phr0, (1-s.phr0)), nrow = 1)

    } else { # put all load in column 1 (active input)
      visib[[1]] <- matrix(c(1,0), nrow = 1)
      phero[[1]] <- matrix(c(1,0), nrow = 1)
    }
    if (ds > 1) {
      for (i in 2:ds) {
        if (i %in% s.indvar) { # distribute load
          visib[[i]] <- do.call(rbind, replicate(ncol(visib[[i-1]]), c(s.vis0, (1-s.vis0)), simplify = F))
          phero[[i]] <- do.call(rbind, replicate(ncol(phero[[i-1]]), c(s.phr0, (1-s.phr0)), simplify = F))

        } else { # put all load in column 1 (active input)
          visib[[i]] <- do.call(rbind, replicate(ncol(visib[[i-1]]), c(1, 0), simplify = F))
          phero[[i]] <- do.call(rbind, replicate(ncol(phero[[i-1]]), c(1, 0), simplify = F))
        }
      }
    }

    # set up initial visibilities and pheromones for functional inputs
    for (i in 1:df) {
      piv <- ds+i
      if (i %in% f.indvar) { # distribute load
        visib[[piv]] <- do.call(rbind, replicate(ncol(visib[[piv-1]]), c(s.vis0, (1-s.vis0)), simplify = F))
        phero[[piv]] <- do.call(rbind, replicate(ncol(phero[[piv-1]]), c(s.phr0, (1-s.phr0)), simplify = F))

      } else { # put all load in column 1 (active input)
        visib[[piv]] <- do.call(rbind, replicate(ncol(visib[[piv-1]]), c(1, 0), simplify = F))
        phero[[piv]] <- do.call(rbind, replicate(ncol(phero[[piv-1]]), c(1, 0), simplify = F))
      }


    }

  } else if(!is.null(fIn)) { # functional-input case ***************************************
  } else if(!is.null(sIn)) { # scalar-input case *******************************************
  } else { # error: no inputs were provided
    stop("The user must provide either a scalar-input matrix, a functional-input list or both of them. None has been detected.")
  }

















  # fill the visibilities and pheromones related to scalar inputs (if there is any)
  if (ds > 0) {
    # the first level only has one possible origin, thus only one row in the matrices
    if (1 %in% s.indvar) {
      visib[[1]] <- matrix(c(s.vis0, (1-s.vis0)), nrow = 1)
      phero[[1]] <- matrix(c(s.phr0, (1-s.phr0)), nrow = 1)
    } else {
      visib[[1]] <- matrix(1)
      phero[[1]] <- matrix(1)
    }

    if (ds > 1) {
      # all other levels might have more than one origin (row), unless the preceding variable is fixed active
      for (i in 2:ds) {
        if (i %in% s.indvar) { # number of columns = 2
          if ((i-1) %in% s.indvar) { # number of rows = 2
            visib[[i]] <- do.call(rbind, replicate(2, c(s.vis0, (1-s.vis0)), simplify = F))
            phero[[i]] <- do.call(rbind, replicate(2, c(s.phr0, (1-s.phr0)), simplify = F))
          } else { # number of rows = 1
            visib[[i]] <- matrix(c(s.phr0, (1-s.vis0)), nrow = 1)
            phero[[i]] <- matrix(c(s.phr0, (1-s.phr0)), nrow = 1)
          }
        } else { # number of columns = 1
          if ((i-1) %in% s.indvar) { # number of rows = 2
            visib[[i]] <- matrix(1, nrow = 2)
            phero[[i]] <- matrix(1, nrow = 2)
          } else { # number of rows = 1
            visib[[i]] <- matrix(1)
            phero[[i]] <- matrix(1)
          }
        }
      }
    }
  }

  # fill the visibilities and pheromones related to functional inputs, one by one (if there is any)
  if (df > 0) {
    # if there is no scalar input in the model, the first level only has one possible origin, thus only one row in the matrices
    if (ds == 0) {
      if (1 %in% f.indvar) { # two columns
        visib[[1]] <- matrix(c(s.vis0, (1-s.vis0)), nrow = 1)
        phero[[1]] <- matrix(c(s.phr0, (1-s.phr0)), nrow = 1)
      } else { # one colum
        visib[[1]] <- matrix(1)
        phero[[1]] <- matrix(1)
      }
    } else {
      if (ds %in% s.indvar) {

      }
      if (1 %in% f.indvar) { # two columns
        visib[[1]] <- matrix(c(s.vis0, (1-s.vis0)), nrow = 1)
        phero[[1]] <- matrix(c(s.phr0, (1-s.phr0)), nrow = 1)
      } else { # one colum
        visib[[1]] <- matrix(1)
        phero[[1]] <- matrix(1)
      }
    }

    if (ds > 1) {
      # all other levels might have more than one origin (row), unless the preceding variable is fixed active
      for (i in 2:ds) {
      }
    }
  }


  # # first fill the initial visibilities and pheromones for the state of scalar inputs (if there is any)
  # if (s.nvar > 0) {
  #   s.acp <- .7 # probability of activating any scalar input
  #   s.ph0 <- .5 # initial pheromones charge for the node where a scalar input is activated
  #
  #   # set the first visibility and pheromones matrices into the list
  #   # those are the only ones with one row because relate the origin with the state of one scalar input from {1,0}
  #   visib[[1]] <- matrix(c(s.acp, (1-s.acp)), nrow = 1)
  #   phero[[1]] <- matrix(c(s.ph0, (1-s.ph0)), nrow = 1)
  #
  #   # set the remaining first visibility and pheromones matrices related to the state of scalar inputs
  #   # those are 2x2 since they go from the one boolean state of the current variable to one boolean state of the next variable
  #   if (s.nvar > 1) {
  #     for (s in 2:s.nvar) {
  #       visib[[s]] <- do.call(rbind, replicate(2, c(s.acp, (1-s.acp)), simplify = F))
  #       phero[[s]] <- do.call(rbind, replicate(2, c(s.ph0,(1-s.ph0)), simplify = F))
  #     }
  #   }
  # }

  # then fill one by one the initial visibilities of the fields related to functional inputs, in thi following order
  # 1. state
  # 2. distance type
  # 3. projection dimension
  # 4. basis family





  # the initial visibilities and pheromones for the state of functional inputs (if there is any)





  print("Well done ;)")



}


checkAndPrint2 <- function(s.state, f.state, f.dim, f.dst, f.fam, k.type) {
  ds <- 2
  df <- 2
  free.fact <- free.values <- fixe.fact <- fixe.values <- c()
  for (i in 1:ds) {
    if (s.state[i] == 0) {
      free.fact <- c(free.fact, paste("State of X", i, sep = ""))
      free.values <- c(free.values, "Inactive, Active")
    } else {
      fixe.fact <- c(fixe.fact, paste("State of X", i, sep = ""))
      fixe.values <- c(fixe.values, "Active")
    }
  }

  for (i in 1:df) {
    if (f.state[i] == 0) {
      free.fact <- c(free.fact, paste("State of F", i, sep = ""))
      free.values <- c(free.values, "Inactive, Active")
    } else {
      fixe.fact <- c(fixe.fact, paste("State of F", i, sep = ""))
      fixe.values <- c(fixe.values, "Active")
    }
  }

  for (i in 1:df) {
    if (length(f.dim[[i]]) > 1) {
      free.fact <- c(free.fact, paste("Dim. of F", i, sep = ""))
      if (0 %in% f.dim[[i]]) {
        free.values <- c(free.values, paste("Orig,", paste(range(f.dim[[i]][-1]), collapse = ":")))
      } else {
        free.values <- c(free.values, paste(range(f.dim[[i]]), collapse = ":"))
      }
    } else {
      fixe.fact <- c(fixe.fact, paste("Dim. of F", i, sep = ""))
      fixe.values <- c(fixe.values, as.character(f.dim[[i]]))
    }
  }

  for (i in 1:df) {
    if (length(f.dst[[i]]) > 1) {
      free.fact <- c(free.fact, paste("Dist. for F", i, sep = ""))
      free.values <- c(free.values, paste(f.dst[[i]], collapse = ", "))
    } else {
      fixe.fact <- c(fixe.fact, paste("Dist. for F", i, sep = ""))
      fixe.values <- c(fixe.values, f.dst[[i]])
    }
  }

  for (i in 1:df) {
    if (length(f.fam[[i]]) > 1) {
      free.fact <- c(free.fact, paste("Basis. for F", i, sep = ""))
      free.values <- c(free.values, paste(f.fam[[i]], collapse = ", "))
    } else {
      fixe.fact <- c(fixe.fact, paste("Basis. for F", i, sep = ""))
      fixe.values <- c(fixe.values, f.fam[[i]])
    }
  }

  if (length(k.type) > 1) {
    free.fact <- c(free.fact, "Kernel type")
    free.values <- c(free.values, paste(k.type, collapse = ", "))
  } else {
    fixe.fact <- c(fixe.fact, "Kernel type")
    fixe.values <- c(fixe.values, k.type)
  }

  cat("Free\n")
  print(knitr::kable(cbind(free.fact, free.values), align = "c", col.names = c("Factor", "Levels"), caption = "Free factors"))
  if (length(fixe.fact) > 0) {
    cat("Fixed\n")
    print(knitr::kable(cbind(fixe.fact, fixe.values), align = "c", col.names = c("Factor", "Levels"), caption = "Fixed factors"))
  }
}
