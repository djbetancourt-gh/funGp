funGpbidon <- function(sIn = NULL, fIn = NULL, sOut, doProj = T, fpDims = NULL, kerType = "matern5_2", disType = "functional",
                       var.hyp = NULL, ls_s.hyp = NULL, ls_f.hyp = NULL, n.starts = 1, n.presample = 20) {
  checkVal_funGp(as.list(environment()))
  # browser()
  # create objects of class funGpProj, funGpKern and funGp
  proj <- new("funGpProj")
  kern <- new("funGpKern")
  model <- new("funGp")

  # extract generic information from user inputs
  sOut <- as.matrix(sOut)
  n.tr <- length(sOut)
  n.presample <- max(n.presample, n.starts)

  # 3 possible cases
  # Case 1: scalar and functional
  # Case 2: functional only
  # Case 3: scalar only
  if (all(!is.null(sIn), !is.null(fIn))) { # Hybrid-input case *******************************************
    # extract information from user inputs specific to the hybrid-input case
    sIn <- as.matrix(sIn)
    ds <- ncol(sIn)
    df <- length(fIn)
    fDims <- sapply(fIn, ncol)

    # Extend to other possible cases!!!!!!!!!!!!!!!!!!
    if (doProj) {
      if (is.null(fpDims)) {
        fpDims <- rep(3, df)
        # fpDims <- c(3,2)
      }

      # project functional inputs
      basis <- fpIn <- J <- list()
      for (i in 1:df) {
        if (fpDims[i] > 0) {
          B <- (eigen(cov(fIn[[i]]))$vectors)[,1:fpDims[i]]
          fpIn[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn[[i]]))
          J[[i]] <- t(B) %*% B
        } else {
          J[[i]] <- B <- diag(ncol(fIn[[i]]))
          fpIn[[i]] <- fIn[[i]]
        }
        basis[[i]] <- B
      }
    } else {
      fpDims <- rep(0, df)
      basis <- J <- lapply(fIn, function(m) diag(ncol(m)))
      fpIn <- fIn
    }

    # compute scalar distance matrices
    sMs <- setScalDistance(sIn, sIn)

    # compute functional distance matrices
    fMs <- setFunDistance(fpIn, fpIn, J)

    # optimize hyperparameters if some if required
    if (!all(is.null(var.hyp), !is.null(ls_s.hyp), !is.null(ls_f.hyp))) {
      varHyp <- var.hyp
      lsHyps <- c(ls_s.hyp, ls_f.hyp)
    } else {
      hypers <- setHypers_SF(sIn, fpIn, J, sMs, fMs, sOut, kerType, n.starts, n.presample)
      varHyp <- hypers[1]
      lsHyps <- hypers[-1]
    }

    # fill funGpKern slots specific to the functional-input case
    kern@s_lsHyps <- lsHyps[1:ds]
    kern@f_lsHyps <- lsHyps[-c(1:ds)]

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_SF(sMs, fMs, sOut, varHyp, lsHyps[1:ds], lsHyps[(ds+1):(ds+df)], kerType)

    # fill funGpProj slots specific to the hybrid-input case
    proj@doProj <- doProj
    proj@fpDims <- fpDims
    proj@basis <- basis
    proj@coefs <- fpIn

    # fill funGp slots specific to the hybrid-input case
    model@ds <- ds
    model@df <- df
    model@fDims <- fDims
    model@sIn <- sIn
    model@fIn <- fIn
    model@type = "hybrid"

  } else if(!is.null(fIn)) { # functional-input case ***************************************
    # extract information from user inputs specific to the functional-input case
    df <- length(fIn)
    fDims <- sapply(fIn, ncol)

    # Extend to other possible cases!!!!!!!!!!!!!!!!!!
    if (all(doProj, is.null(fpDims))) {
      fpDims <- rep(3, df)
      # fpDims <- c(3,2)
    }

    # project functional inputs
    basis <- fpIn <- J <- list()
    for (i in 1:df) { # functional-input case
      B <- (eigen(cov(fIn[[i]]))$vectors)[,1:fpDims[i]]
      fpIn[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn[[i]]))
      J[[i]] <- t(B) %*% B
      basis[[i]] <- B
    }

    # compute functional distance matrices
    fMs <- setFunDistance(fpIn, fpIn, J)

    # optimize hyperparameters if some if required
    if (!all(is.null(var.hyp), !is.null(ls_f.hyp))) {
      varHyp <- var.hyp
      lsHyps <- ls_f.hyp
    } else {
      hypers <- setHypers_F(fpIn, J, fMs, sOut, kerType, n.starts, n.presample)
      varHyp <- hypers[1]
      lsHyps <- hypers[-1]
    }

    # fill funGpKern slots specific to the functional-input case
    kern@f_lsHyps <- lsHyps

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_F(fMs, sOut, varHyp, lsHyps, kerType)

    # fill funGpProj slots specific to the hybrid-input case
    proj@doProj <- doProj
    proj@fpDims <- fpDims
    proj@basis <- basis
    proj@coefs <- fpIn

    # fill funGp slots specific to the functional-input case
    model@ds <- 0
    model@df <- df
    model@fDims <- fDims
    model@fIn <- fIn
    model@type = "functional"

  } else if(!is.null(sIn)) { # scalar-input case *******************************************
    # extract information from user inputs specific to the scalar-input case
    sIn <- as.matrix(sIn)
    ds <- ncol(sIn)

    # compute scalar distance matrices
    sMs <- setScalDistance(sIn, sIn)

    # optimize hyperparameters if some if required
    if (!all(is.null(var.hyp), !is.null(ls_f.hyp))) {
      varHyp <- var.hyp
      lsHyps <- ls_s.hyp
    } else {
      hypers <- setHypers_S(sIn, sMs, sOut, kerType, n.starts, n.presample)
      varHyp <- hypers[1]
      lsHyps <- hypers[-1]
    }

    # fill funGpKern slots specific to the scalar-input case
    kern@s_lsHyps <- lsHyps

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_S(sMs, sOut, varHyp, lsHyps, kerType)

    # fill funGp slots specific to the scalar-input case
    model@ds <- ds
    model@df <- 0
    model@sIn <- sIn
    model@type = "scalar"

  } else { # error: no inputs were provided
    stop("User must provide either a scalar-input matrix, a functional-input list or both of them. None has been detected.")
  }

  # fill general funGpKern slots
  kern@kerType <- kerType
  kern@disType <- disType
  kern@varHyp <- varHyp

  # fill general funGpModel slots
  model@call <- match.call()
  model@sOut <- sOut
  model@n.tot <- n.tr
  model@n.tr <- n.tr
  model@proj <- proj
  model@kern <- kern

  return(model)
}
# -------------------------------------------------------------------------------------------------------------------------------------


# ==========================================================================================================
# Updating functions
# ==========================================================================================================

# Function to delete some data
# ----------------------------------------------------------------------------------------------------------
upd_del <- function(model, ind.dl, remake = F) {
  # duplicate the original model to build the updated one
  modelup <- model

  # check for validty of substituting data
  ind.dl <- check_del(as.list(environment()))

  if (model@type == "hybrid") { # hybrid-input case *******************************************
    # extract inputs from original model and remove points according to deletion indices
    sIn <- model@sIn[-ind.dl,,drop = F]
    fIn <- lapply(model@fIn, function(M) M[-ind.dl,])
    sOut <- model@sOut[-ind.dl,,drop = F]

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      modelup <- funGp(sIn = sIn, fIn = fIn, sOut = sOut, doProj = model@proj@doProj, fpDims = model@proj@fpDims,
                       kerType = model@kern@kerType, disType = model@kern@disType, var.hyp = model@kern@varHyp,
                       ls_s.hyp = model@kern@s_lsHyps, ls_f.hyp = model@kern@f_lsHyps)
    } else {
      modelup@sIn <- sIn
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }

  } else if (model@type == "functional") { # functional-input case *******************************************
    # extract inputs from original model
    fIn <- model@fIn

    # remove points according to deletion indices
    fIn <- lapply(fIn, function(M) M[-ind.dl,])
    sOut <- sOut[-ind.dl,,drop = F]

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      modelup <- funGp(fIn = fIn, sOut = sOut, doProj = model@proj@doProj, fpDims = model@proj@fpDims,
                       kerType = model@kern@kerType, disType = model@kern@disType,
                       var.hyp = model@kern@varHyp, ls_f.hyp = model@kern@f_lsHyps)
    } else {
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }

  } else { # scalar-input case *******************************************
    # extract inputs from original model
    sIn <- model@sIn

    # remove points according to deletion indices
    sIn <- sIn[-ind.dl,,drop = F]
    sOut <- sOut[-ind.dl,,drop = F]

    # request new model to refunGp if requested
    if (remake) {
      modelup <- funGp(sIn = sIn, sOut = sOut, kerType = model@kern@kerType,
                       var.hyp = model@kern@varHyp, ls_s.hyp = model@kern@s_lsHyps)
    } else {
      modelup@sIn <- sIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }
  }

  return(modelup)
}
# ----------------------------------------------------------------------------------------------------------


# Function to substitute some data
# ----------------------------------------------------------------------------------------------------------
upd_subData <- function(model, sIn.sb, fIn.sb, sOut.sb, ind.sb, remake = F) {
  # duplicate the original model to build the updated one
  modelup <- model

  # extract generic information from the model
  sOut <- model@sOut

  # identify the special case of only substituting in sOut
  if(all(is.null(sIn.sb), is.null(fIn.sb), !is.null(sOut.sb))) justOut <- T else justOut <- F

  # provide substituting output if not specified by the user
  if(is.null(sOut.sb)) sOut.sb <- sOut[ind.sb,,drop = F]

  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # extract inputs from original model
    sIn <- model@sIn
    fIn <- model@fIn

    # provide substituting inputs if not specified by the user
    if(is.null(sIn.sb)) sIn.sb <- sIn[ind.sb,,drop = F]
    if(is.null(fIn.sb)) fIn.sb <- lapply(fIn, function(M) M[ind.sb,,drop = F])

    # check for validty of substituting data
    check_sub(as.list(environment()))

    # check for duplicates in the substituting points
    ind.dp <- check_duplicates_SF(sIn.sb, fIn.sb, sIn.sb, fIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
                    "Substitution is skipped.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.sb <- sIn.sb[-ind.dp,,drop = F]
      fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = F])
      sOut.sb <- sOut.sb[-ind.dp,,drop = F]
      ind.sb <- ind.sb[-ind.dp]
    }

    # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
    sIn.exsb <- sIn[-ind.sb,]
    fIn.exsb <- lapply(fIn, function(M) M[-ind.sb,,drop = F])
    ind.dp <- check_duplicates_SF(sIn.exsb, fIn.exsb, sIn.sb, fIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
                    "contained in the model. The model is returned in its original state.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.sb <- sIn.sb[-ind.dp,,drop = F]
      fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = F])
      Sout.sb <- sOut.sb[-ind.dp,,drop = F]
      ind.sb <- ind.sb[-ind.dp]
    }

    # recover inputs and outputs after duplicates check
    sIn[ind.sb,] <- sIn.sb
    fIn <- mapply(function(M, x) {M[ind.sb,] <- x; return(M)}, fIn, fIn.sb)
    sOut[ind.sb,] <- sOut.sb

    # request new model to refunGp if requested
    if (remake) {
      if (justOut) {
        modelup@preMats$LInvY <- backsolve(model@preMats$L, sOut, upper.tri = F)
        modelup@sIn <- sIn
        modelup@fIn <- fIn
        modelup@sOut <- sOut
        model@n.tot <- length(sOut)
      } else {
        modelup <- funGp(sIn = sIn, fIn = fIn, sOut = sOut, doProj = model@proj@doProj, fpDims = model@proj@fpDims,
                         kerType = model@kern@kerType, disType = model@kern@disType, var.hyp = model@kern@varHyp,
                         ls_s.hyp = model@kern@s_lsHyps, ls_f.hyp = model@kern@f_lsHyps)
      }
    } else {
      modelup@sIn <- sIn
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }

  } else if (model@type == "functional") { # functional-input case *******************************************
    # extract inputs from original model
    fIn <- model@fIn

    # provide substituting inputs if not specified by the user
    if(is.null(fIn.sb)) fIn.sb <- lapply(fIn, function(M) M[ind.sb,,drop = F])

    # check for validty of substituting data
    check_sub(as.list(environment()))

    # check for duplicates in the substituting points
    ind.dp <- check_duplicates_F(fIn.sb, fIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
                    "Substitution is skipped.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = F])
      sOut.sb <- sOut.sb[-ind.dp,,drop = F]
      ind.sb <- ind.sb[-ind.dp]
    }

    # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
    fIn.exsb <- lapply(fIn, function(M) M[-ind.sb,,drop = F])
    ind.dp <- check_duplicates_F(fIn.exsb, fIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
                    "contained in the model. The model is returned in its original state.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = F])
      Sout.sb <- sOut.sb[-ind.dp,,drop = F]
      ind.sb <- ind.sb[-ind.dp]
    }

    # recover inputs and outputs after duplicates check
    fIn <- mapply(function(M, x) {M[ind.sb,] <- x; return(M)}, fIn, fIn.sb)
    sOut[ind.sb,] <- sOut.sb

    # request new model to refunGp if requested
    if (remake) {
      if (justOut) {
        modelup@preMats$LInvY <- backsolve(model@preMats$L, sOut, upper.tri = F)
        modelup@fIn <- fIn
        modelup@sOut <- sOut
        model@n.tot <- length(sOut)
      } else {
        modelup <- funGp(fIn = fIn, sOut = sOut, doProj = model@proj@doProj, fpDims = model@proj@fpDims,
                         kerType = model@kern@kerType, disType = model@kern@disType,
                         var.hyp = model@kern@varHyp, ls_f.hyp = model@kern@f_lsHyps)
      }
    } else {
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }

  } else { # scalar-input case *******************************************
    # extract inputs from original model
    sIn <- model@sIn

    # provide substituting inputs if not specified by the user
    if(is.null(sIn.sb)) sIn.sb <- sIn[ind.sb]

    # check for validty of substituting data
    check_sub(as.list(environment()))

    # check for duplicates in the substituting points
    ind.dp <- check_duplicates_S(sIn.sb, sIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
                    "Substitution is skipped.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.sb <- sIn.sb[-ind.dp,,drop = F]
      sOut.sb <- sOut.sb[-ind.dp,,drop = F]
      ind.sb <- ind.sb[-ind.dp]
    }

    # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
    sIn.exsb <- sIn[-ind.sb,,drop = F]
    ind.dp <- check_duplicates_S(sIn.exsb, sIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
                    "contained in the model. The model is returned in its original state.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.sb <- sIn.sb[-ind.dp,,drop = F]
      sOut.sb <- sOut.sb[-ind.dp,,drop = F]
      ind.sb <- ind.sb[-ind.dp]
    }

    # recover inputs and outputs after duplicates check
    sIn[ind.sb,] <- sIn.sb
    sOut[ind.sb,] <- sOut.sb

    # request new model to refunGp if requested
    if (remake) {
      if (justOut) {
        modelup@preMats$LInvY <- backsolve(model@preMats$L, sOut, upper.tri = F)
        modelup@sIn <- sIn
        modelup@sOut <- sOut
        model@n.tot <- length(sOut)
      } else {
        modelup <- funGp(sIn = sIn, sOut = sOut, kerType = model@kern@kerType,
                         var.hyp = model@kern@varHyp, ls_s.hyp = model@kern@s_lsHyps)
      }
    } else {
      modelup@sIn <- sIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }
  }

  return(modelup)
}
# ----------------------------------------------------------------------------------------------------------


# Function to add some data
# ----------------------------------------------------------------------------------------------------------
upd_add <- function(model, sIn.nw, fIn.nw, sOut.nw, remake = F) {
  # browser()
  # check validty of substituting data
  check_add(as.list(environment()))

  # duplicate the original model to build the updated one
  modelup <- model

  # extract generic information from the model
  sOut <- model@sOut

  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # extract inputs from original model
    sIn <- model@sIn
    fIn <- model@fIn

    # check for duplicates in the new points
    ind.dp <- check_duplicates_SF(sIn.nw, fIn.nw, sIn.nw, fIn.nw)

    if (length(sOut) == length(ind.dp)) {
      warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
                    "Substitution is skipped.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.nw <- sIn.nw[-ind.dp,,drop = F]
      fIn.nw <- lapply(fIn.nw, function(M) M[-ind.dp,,drop = F])
      sOut.nw <- sOut.nw[-ind.dp,,drop = F]
    }

    # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
    ind.dp <- check_duplicates_SF(sIn, fIn, sIn.nw, fIn.nw)

    if (length(sOut) == length(ind.dp)) {
      warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
                    "contained in the model. The model is returned in its original state.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.nw <- sIn.nw[-ind.dp,,drop = F]
      fIn.nw <- lapply(fIn.nw, function(M) M[-ind.dp,,drop = F])
      sOut.nw <- sOut.nw[-ind.dp,,drop = F]
    }

    # recover inputs and outputs after duplicates check
    sIn <- rbind(sIn, sIn.nw)
    fIn <- Map(rbind, fIn, fIn.nw)
    sOut <- rbind(sOut, sOut.nw)

    # request new model to refunGp if requested
    if (remake) {
      modelup <- funGpbidon(sIn = sIn, fIn = fIn, sOut = sOut, var.hyp = model@kern@varHyp,
                            ls_s.hyp = model@kern@s_lsHyps, ls_f.hyp = model@kern@f_lsHyps)
    } else {
      modelup@sIn <- sIn
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }

  } else if (model@type == "functional") { # functional-input case *******************************************
    # extract inputs from original model
    fIn <- model@fIn

    # check for duplicates in the new points
    ind.dp <- check_duplicates_F(fIn.nw, fIn.nw)

    if (length(sOut) == length(ind.dp)) {
      warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
                    "Substitution is skipped.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      fIn.nw <- lapply(fIn.nw, function(M) M[-ind.dp,,drop = F])
      sOut.nw <- sOut.nw[-ind.dp,,drop = F]
    }

    # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
    ind.dp <- check_duplicates_F(fIn, fIn.nw)

    if (length(sOut) == length(ind.dp)) {
      warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
                    "contained in the model. The model is returned in its original state.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      fIn.nw <- lapply(fIn.nw, function(M) M[-ind.dp,,drop = F])
      sOut.nw <- sOut.nw[-ind.dp,,drop = F]
    }

    # recover inputs and outputs after duplicates check
    fIn <- Map(rbind, fIn, fIn.nw)
    sOut <- rbind(sOut, sOut.nw)

    # request new model to refunGp if requested
    if (remake) {
      modelup <- funGpbidon(fIn = fIn, sOut = sOut, var.hyp = model@kern@varHyp, ls_f.hyp = model@kern@f_lsHyps)
    } else {
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }

  } else { # scalar-input case *******************************************
    # extract inputs from original model
    sIn <- model@sIn

    # check for duplicates in the new points
    ind.dp <- check_duplicates_S(sIn.nw, sIn.nw)

    if (length(sOut) == length(ind.dp)) {
      warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
                    "Substitution is skipped.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.nw <- sIn.nw[-ind.dp,,drop = F]
      sOut.nw <- sOut.nw[-ind.dp,,drop = F]
    }

    # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
    ind.dp <- check_duplicates_S(sIn, sIn.nw)

    if (length(sOut) == length(ind.dp)) {
      warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
                    "contained in the model. The model is returned in its original state.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.nw <- sIn.nw[-ind.dp,,drop = F]
      sOut.nw <- sOut.nw[-ind.dp,,drop = F]
    }

    # recover inputs and outputs after duplicates check
    sIn <- rbind(sIn, sIn.nw)
    sOut <- rbind(sOut, sOut.nw)

    # request new model to refunGp if requested
    if (remake) {
      modelup <- funGpbidon(sIn = sIn, sOut = sOut, var.hyp = model@kern@varHyp, ls_s.hyp = model@kern@s_lsHyps)
    } else {
      modelup@sIn <- sIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }
  }

  return(modelup)
}
# ----------------------------------------------------------------------------------------------------------


# Function to substitute some Hyperparameters
# ----------------------------------------------------------------------------------------------------------
upd_subHypers <- function(model, var.sb, ls_s.sb, ls_f.sb) {
  # browser()
  # check validty of substituting hypers
  check_add(as.list(environment()))

  # duplicate the original model to build the updated one
  modelup <- model

  if (model@type == "hybrid") { # Hybrid-input case *******************************************


    # remake model according to the case
    # 1. var, ls_s, ls_f
    # 2. var, ls_s
    # 3. var, ls_f
    # 4. ls_s, ls_f
    # 5. ls_s
    # 6. ls_f
    if (!is.null(var.sb)) {
      if (!is.null(ls_s.sb)) {
        if (!is.null(ls_f.sb)) { # case 1: Substitute all hypers

        } else { # case 2: Substitute the var and ls_s

        }
      } else if (!is.null(ls_f.sb)) { # case 3: Substitute the var and ls_f

      }
    } else if (!is.null(ls_s.sb)) {
      if (!is.null(ls_f.sb)) { # case 4: Substitute ls_s and ls_f

      } else { # case 5: Substitute only ls_s

      }
    } else { # case 6: Substitute only ls_f

    }


    # if (all(!is.null(var.sb), !is.null(ls_s.sb), !is.null(ls_f.sb))) {
    #   modelup <- upd_add(model = object, var.sb = var.sb, ls_s.sb = ls_s.sb, ls_f.sb = ls_f.sb)
    # } else if (all(!is.null(var.sb), !is.null(ls_s.sb))) {
    #   modelup <- upd_add(model = object, var.sb = var.sb, ls_s.sb = ls_s.sb)
    # } else if (all(!is.null(var.sb), !is.null(ls_f.sb))) {
    #   modelup <- upd_add(model = object, var.sb = var.sb, ls_f.sb = ls_f.sb)
    # } else if (all(!is.null(ls_s.sb), !is.null(ls_f.sb))) {
    #   modelup <- upd_add(model = object, ls_s.sb = ls_s.sb, ls_f.sb = ls_f.sb)
    # } else if (!is.null(var.sb)) {
    #   modelup <- upd_add(model = object, var.sb = var.sb)
    # } else if (!is.null(ls_s.sb)) {
    #   modelup <- upd_add(model = object, ls_s.sb = ls_s.sb)
    # } else if (!is.null(ls_f.sb)) {
    #   modelup <- upd_add(model = object, ls_f.sb = ls_f.sb)
    # }

  } else if (model@type == "functional") { # functional-input case *******************************************


  } else { # scalar-input case *******************************************

  }

  return(modelup)








}
# ----------------------------------------------------------------------------------------------------------

# ==========================================================================================================
# Validators of updating functions
# ==========================================================================================================

# Function to check if the parameters for data deletion are ok
# ----------------------------------------------------------------------------------------------------------
check_del <- function(env) {
  # recover the model
  model <- env$model

  # is there any duplicate?
  if (any(duplicated(env$ind.dl))) {
    warning("Some elements in the deletion index vector ind.dl are duplicated. Duplicates are dropped.")
    env$ind.dl <- unique(env$ind.dl)
  }

  # are indices in the correct range?
  if (any(env$ind.dl < 1) || any(env$ind.dl > model@n.tot)) {
    stop("Some elements in the detelion index vector ind.dl are not in [1, model@n.tot]. Please check your vector.")
  }

  # is the maximum number of points to delete respected?
  if (model@n.tot - length(env$ind.dl) < 2) {
    stop("The number of points to delete cannot exceed model@n.tot - 2. Please check your vector.")
  }
  return(env$ind.dl)
}
# ----------------------------------------------------------------------------------------------------------


# Function to check if the parameters for data substitution are ok
# ----------------------------------------------------------------------------------------------------------
check_sub <- function(env) {
  # browser()
  # recover the model
  model <- env$model
  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # consistency in number of points
    if (!check_nrows_mlm(env$sIn.sb, env$fIn.sb, env$sOut.sb)) {
      stop("Inconsistent number of points. Please check that sIn.sb, each matrix in fIn.sb and sOut.sb have all the same number\nof rows.")
    }
    if (nrow(env$sIn.sb) != length(env$ind.sb)) {
      stop(paste("Inconsistent number of points in your replacement index vector ind.sb. Please check that it has the same number of\n",
                 "rows than sIn.sb, each matrix in fIn.sb and sOut.sb.", sep = ""))
    }

  } else if (model@type == "functional") {
    # consistency in number of points
    if (!check_nrows_lm(env$fIn.sb, env$sOut.sb)) {
      stop("Inconsistent number of points. Please check that each matrix in fIn.sb and sOut.sb have all the same number of rows.")
    }
    if (nrow(env$fIn.sb[[1]]) != length(env$ind.sb)) {
      stop(paste("Inconsistent number of points in your replacement index vector ind.sb. Please check that it has the same number of\n",
                 "rows than each matrix in fIn.sb and sOut.sb.", sep = ""))
    }

  } else if (model@type == "scalar") {
    # consistency in number of points
    if (!check_nrows_mm(env$sIn.sb, env$sOut.sb)) {
      stop("Inconsistent number of points. Please check that sIn.sb and sOut.sb have the same number of rows.")
    }
    if (nrow(env$sIn.sb) != length(env$ind.sb)) {
      stop(paste("Inconsistent number of points in your replacement index vector ind.sb. Please check that it has the same number of\n",
                 "rows than sIn.sb and sOut.sb.", sep = ""))
    }
  }

  # validity of subtituting index
  # should be an integer
  if (any(env$ind.sb %% 1 != 0)) {
    stop(paste(c("Replacement indices should be integer numbers. Please check the following positions of your ind.sb vector: ",
                 which(env$ind.sb %% 1 != 0)), sep = "", collapse = " "))
  }
  # should be in [1, n.tot]
  if (any(env$ind.dl < 1) || any(env$ind.dl > model@n.tot)) {
    stop(paste(c("Replacement indices shuold be integers in [1, model@n.tot]. Please check the following positions of your ind.sb vector: ",
                 which(!(env$ind.sb > 0 & env$ind.sb < model@n.tr))), sep = "", collapse = " "))
  }
  # should not contain duplicates
  if (anyDuplicated(env$ind.sb)) {
    stop(paste(c("Replacement indices shuold be unique. Please check the following positions of your ind.sb vector: ",
                 which(duplicated(env$ind.sb) | duplicated(env$ind.sb, fromLast = TRUE))), sep = "", collapse = " "))
  }
}
# ----------------------------------------------------------------------------------------------------------


# Function to check if the parameters for data addition are ok
# ----------------------------------------------------------------------------------------------------------
check_add <- function(env) {
  # browser()
  # recover the model
  model <- env$model
  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # consistency in number of points
    if (!check_nrows_mlm(env$sIn.nw, env$fIn.nw, env$sOut.nw)) {
      stop("Inconsistent number of points. Please check that sIn.nw, each matrix in fIn.nw and sOut.nw have all the same number\nof rows.")
    }

  } else if (model@type == "functional") {
    # consistency in number of points
    if (!check_nrows_lm(env$fIn.nw, env$sOut.nw)) {
      stop("Inconsistent number of points. Please check that each matrix in fIn.nw and sOut.nw have all the same number of rows.")
    }

  } else if (model@type == "scalar") {
    # consistency in number of points
    if (!check_nrows_mm(env$sIn.nw, env$sOut.nw)) {
      stop("Inconsistent number of points. Please check that sIn.nw and sOut.nw have the same number of rows.")
    }
  }
}
# ----------------------------------------------------------------------------------------------------------
