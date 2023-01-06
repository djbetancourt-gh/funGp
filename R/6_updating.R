# ==========================================================================================================
# Updating functions
# ==========================================================================================================

# Function to delete some data
# ----------------------------------------------------------------------------------------------------------
upd_del <- function(model, ind.dl, remake, trace, pbars, control.optim) {
  # check for validty of substituting data
  ind.dl <- check_del(as.list(environment()))

  if (model@type == "hybrid") { # hybrid-input case *******************************************
    # extract inputs from original model and remove points according to deletion indices
    sIn <- model@sIn[-ind.dl,,drop = FALSE]
    fIn <- lapply(model@fIn, function(M) M[-ind.dl,])
    sOut <- model@sOut[-ind.dl,,drop = FALSE]

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      modelup <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut, kerType = model@kern@kerType,
                      f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims,
                      f_basType = model@f_proj@basType, var.hyp = model@kern@varHyp,
                      ls_s.hyp = model@kern@s_lsHyps, ls_f.hyp = model@kern@f_lsHyps,
                      trace = trace, pbars = pbars, control.optim = control.optim)
    } else {
      modelup <- model
      # update projection
      bcj <- dimReduction(fIn, model@df, model@proj@pdims, model@f_proj@basType)
      modelup@f_proj@basis <- bcj$basis
      modelup@f_proj@coefs <- bcj$coefs
      # update data points
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
    sOut <- sOut[-ind.dl,,drop = FALSE]

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      modelup <- fgpm(fIn = fIn, sOut = sOut,
                      kerType = model@kern@kerType, f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims,
                      f_basType = model@f_proj@basType, var.hyp = model@kern@varHyp, ls_f.hyp = model@kern@f_lsHyps,
                      trace = trace, pbars = pbars, control.optim = control.optim)
    } else {
      modelup <- model
      # update projection
      bcj <- dimReduction(fIn, model@df, model@proj@pdims, model@f_proj@basType)
      modelup@f_proj@basis <- bcj$basis
      modelup@f_proj@coefs <- bcj$coefs
      # update data points
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }

  } else { # scalar-input case *******************************************
    # extract inputs from original model
    sIn <- model@sIn

    # remove points according to deletion indices
    sIn <- sIn[-ind.dl,,drop = FALSE]
    sOut <- sOut[-ind.dl,,drop = FALSE]

    # request new model to fgpm if indicated
    if (remake) {
      modelup <- fgpm(sIn = sIn, sOut = sOut, kerType = model@kern@kerType,
                      var.hyp = model@kern@varHyp, ls_s.hyp = model@kern@s_lsHyps,
                      trace = trace, pbars = pbars, control.optim = control.optim)
    } else {
      modelup <- model
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
upd_subData <- function(model, sIn.sb, fIn.sb, sOut.sb, ind.sb, remake, trace, pbars, control.optim) {
  # extract generic information from the model
  sOut <- model@sOut

  # identify the special case of only substituting in sOut
  if(all(is.null(sIn.sb), is.null(fIn.sb), !is.null(sOut.sb))) justOut <- T else justOut <- F

  # provide substituting output if not specified by the user
  if(is.null(sOut.sb)) sOut.sb <- sOut[ind.sb,,drop = FALSE]

  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # extract inputs from original model
    sIn <- model@sIn
    fIn <- model@fIn

    # provide substituting inputs if not specified by the user
    if(is.null(sIn.sb)) sIn.sb <- sIn[ind.sb,,drop = FALSE]
    if(is.null(fIn.sb)) fIn.sb <- lapply(fIn, function(M) M[ind.sb,,drop = FALSE])

    # check for validty of substituting data
    check_subData(as.list(environment()))

    # check for duplicates in the substituting points
    ind.dp <- check_duplicates_SF(sIn.sb, fIn.sb, sIn.sb, fIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
                    "Substitution is skipped.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.sb <- sIn.sb[-ind.dp,,drop = FALSE]
      fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = FALSE])
      sOut.sb <- sOut.sb[-ind.dp,,drop = FALSE]
      ind.sb <- ind.sb[-ind.dp]
    }

    # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
    sIn.exsb <- sIn[-ind.sb,]
    fIn.exsb <- lapply(fIn, function(M) M[-ind.sb,,drop = FALSE])
    ind.dp <- check_duplicates_SF(sIn.exsb, fIn.exsb, sIn.sb, fIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
                    "contained in the model. The model is returned in its original state.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.sb <- sIn.sb[-ind.dp,,drop = FALSE]
      fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = FALSE])
      Sout.sb <- sOut.sb[-ind.dp,,drop = FALSE]
      ind.sb <- ind.sb[-ind.dp]
    }

    # recover inputs and outputs after duplicates check
    sIn[ind.sb,] <- sIn.sb
    fIn <- mapply(function(M, x) {M[ind.sb,] <- x; return(M)}, fIn, fIn.sb)
    sOut[ind.sb,] <- sOut.sb

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      if (justOut) {
        modelup <- model
        modelup@preMats$LInvY <- backsolve(model@preMats$L, sOut, upper.tri = FALSE)
        modelup@sIn <- sIn
        modelup@fIn <- fIn
        modelup@sOut <- sOut
        model@n.tot <- length(sOut)
      } else {
        modelup <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut, kerType = model@kern@kerType,
                        f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims,
                        f_basType = model@f_proj@basType, var.hyp = model@kern@varHyp,
                        ls_s.hyp = model@kern@s_lsHyps, ls_f.hyp = model@kern@f_lsHyps,
                        trace = trace, pbars = pbars, control.optim = control.optim)
      }
    } else {
      modelup <- model
      # update projection
      bcj <- dimReduction(fIn, model@df, model@proj@pdims, model@f_proj@basType)
      modelup@f_proj@basis <- bcj$basis
      modelup@f_proj@coefs <- bcj$coefs
      # update data points
      modelup@sIn <- sIn
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      model@n.tot <- length(sOut)
    }

  } else if (model@type == "functional") { # functional-input case *******************************************
    # extract inputs from original model
    fIn <- model@fIn

    # provide substituting inputs if not specified by the user
    if(is.null(fIn.sb)) fIn.sb <- lapply(fIn, function(M) M[ind.sb,,drop = FALSE])

    # check for validty of substituting data
    check_subData(as.list(environment()))

    # check for duplicates in the substituting points
    ind.dp <- check_duplicates_F(fIn.sb, fIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
                    "Substitution is skipped.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = FALSE])
      sOut.sb <- sOut.sb[-ind.dp,,drop = FALSE]
      ind.sb <- ind.sb[-ind.dp]
    }

    # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
    fIn.exsb <- lapply(fIn, function(M) M[-ind.sb,,drop = FALSE])
    ind.dp <- check_duplicates_F(fIn.exsb, fIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
                    "contained in the model. The model is returned in its original state.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = FALSE])
      Sout.sb <- sOut.sb[-ind.dp,,drop = FALSE]
      ind.sb <- ind.sb[-ind.dp]
    }

    # recover inputs and outputs after duplicates check
    fIn <- mapply(function(M, x) {M[ind.sb,] <- x; return(M)}, fIn, fIn.sb)
    sOut[ind.sb,] <- sOut.sb

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      if (justOut) {
        modelup <- model
        modelup@preMats$LInvY <- backsolve(model@preMats$L, sOut, upper.tri = FALSE)
        modelup@fIn <- fIn
        modelup@sOut <- sOut
        model@n.tot <- length(sOut)
      } else {
        modelup <- fgpm(fIn = fIn, sOut = sOut,
                        kerType = model@kern@kerType, f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims,
                        f_basType = model@f_proj@basType, var.hyp = model@kern@varHyp, ls_f.hyp = model@kern@f_lsHyps,
                        trace = trace, pbars = pbars, control.optim = control.optim)
      }
    } else {
      modelup <- model
      # update projection
      bcj <- dimReduction(fIn, model@df, model@proj@pdims, model@f_proj@basType)
      modelup@f_proj@basis <- bcj$basis
      modelup@f_proj@coefs <- bcj$coefs
      # update data points
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
    check_subData(as.list(environment()))

    # check for duplicates in the substituting points
    ind.dp <- check_duplicates_S(sIn.sb, sIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
                    "Substitution is skipped.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.sb <- sIn.sb[-ind.dp,,drop = FALSE]
      sOut.sb <- sOut.sb[-ind.dp,,drop = FALSE]
      ind.sb <- ind.sb[-ind.dp]
    }

    # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
    sIn.exsb <- sIn[-ind.sb,,drop = FALSE]
    ind.dp <- check_duplicates_S(sIn.exsb, sIn.sb)

    if (length(ind.sb) == length(ind.dp)) {
      warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
                    "contained in the model. The model is returned in its original state.", sep = ""))
      return(model)
    } else if (length(ind.dp) > 0) {
      warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
                    "Duplicate substitute points: ", ind.dp, sep = ""))
      sIn.sb <- sIn.sb[-ind.dp,,drop = FALSE]
      sOut.sb <- sOut.sb[-ind.dp,,drop = FALSE]
      ind.sb <- ind.sb[-ind.dp]
    }

    # recover inputs and outputs after duplicates check
    sIn[ind.sb,] <- sIn.sb
    sOut[ind.sb,] <- sOut.sb

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      if (justOut) {
        modelup <- model
        modelup@preMats$LInvY <- backsolve(model@preMats$L, sOut, upper.tri = FALSE)
        modelup@sIn <- sIn
        modelup@sOut <- sOut
        model@n.tot <- length(sOut)
      } else {
        modelup <- fgpm(sIn = sIn, sOut = sOut, kerType = model@kern@kerType,
                        var.hyp = model@kern@varHyp, ls_s.hyp = model@kern@s_lsHyps,
                        trace = trace, pbars = pbars, control.optim = control.optim)
      }
    } else {
      modelup <- model
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
upd_add <- function(model, sIn.nw, fIn.nw, sOut.nw, remake, trace, pbars, control.optim) {
  # check validty of substituting data
  check_add(as.list(environment()))

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
      sIn.nw <- sIn.nw[-ind.dp,,drop = FALSE]
      fIn.nw <- lapply(fIn.nw, function(M) M[-ind.dp,,drop = FALSE])
      sOut.nw <- sOut.nw[-ind.dp,,drop = FALSE]
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
      sIn.nw <- sIn.nw[-ind.dp,,drop = FALSE]
      fIn.nw <- lapply(fIn.nw, function(M) M[-ind.dp,,drop = FALSE])
      sOut.nw <- sOut.nw[-ind.dp,,drop = FALSE]
    }

    # recover inputs and outputs after duplicates check
    sIn <- rbind(sIn, sIn.nw)
    fIn <- Map(rbind, fIn, fIn.nw)
    sOut <- rbind(sOut, sOut.nw)

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      modelup <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut, kerType = model@kern@kerType,
                      f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims,
                      f_basType = model@f_proj@basType, var.hyp = model@kern@varHyp,
                      ls_s.hyp = model@kern@s_lsHyps, ls_f.hyp = model@kern@f_lsHyps,
                      trace = trace, pbars = pbars, control.optim = control.optim)
    } else {
      modelup <- model
      # update projection
      bcj <- dimReduction(fIn, model@df, model@proj@pdims, model@f_proj@basType)
      modelup@f_proj@basis <- bcj$basis
      modelup@f_proj@coefs <- bcj$coefs
      # update data points
      modelup@sIn <- sIn
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      modelup@n.tot <- length(sOut)
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
      fIn.nw <- lapply(fIn.nw, function(M) M[-ind.dp,,drop = FALSE])
      sOut.nw <- sOut.nw[-ind.dp,,drop = FALSE]
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
      fIn.nw <- lapply(fIn.nw, function(M) M[-ind.dp,,drop = FALSE])
      sOut.nw <- sOut.nw[-ind.dp,,drop = FALSE]
    }

    # recover inputs and outputs after duplicates check
    fIn <- Map(rbind, fIn, fIn.nw)
    sOut <- rbind(sOut, sOut.nw)

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      modelup <- fgpm(fIn = fIn, sOut = sOut,
                      kerType = model@kern@kerType, f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims,
                      f_basType = model@f_proj@basType, var.hyp = model@kern@varHyp, ls_f.hyp = model@kern@f_lsHyps,
                      trace = trace, pbars = pbars, control.optim = control.optim)
    } else {
      modelup <- model
      # update projection
      bcj <- dimReduction(fIn, model@df, model@proj@pdims, model@f_proj@basType)
      modelup@f_proj@basis <- bcj$basis
      modelup@f_proj@coefs <- bcj$coefs
      # update data points
      modelup@fIn <- fIn
      modelup@sOut <- sOut
      modelup@n.tot <- length(sOut)
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
      sIn.nw <- sIn.nw[-ind.dp,,drop = FALSE]
      sOut.nw <- sOut.nw[-ind.dp,,drop = FALSE]
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
      sIn.nw <- sIn.nw[-ind.dp,,drop = FALSE]
      sOut.nw <- sOut.nw[-ind.dp,,drop = FALSE]
    }

    # recover inputs and outputs after duplicates check
    sIn <- rbind(sIn, sIn.nw)
    sOut <- rbind(sOut, sOut.nw)

    # the model is re-made if this is the last one in the sequence of requested tasks
    if (remake) {
      modelup <- fgpm(sIn = sIn, sOut = sOut, kerType = model@kern@kerType,
                      var.hyp = model@kern@varHyp, ls_s.hyp = model@kern@s_lsHyps,
                      trace = trace, pbars = pbars, control.optim = control.optim)
    } else {
      modelup <- model
      modelup@sIn <- sIn
      modelup@sOut <- sOut
      modelup@n.tot <- length(sOut)
    }
  }

  return(modelup)
}
# ----------------------------------------------------------------------------------------------------------


# Function to substitute some Hyperparameters
# ----------------------------------------------------------------------------------------------------------
upd_subHypers <- function(model, var.sb, ls_s.sb, ls_f.sb, trace, pbars, control.optim) {
  # check validty of substituting hypers
  check_subHypers(as.list(environment()))

  # var is always necessary, so if no specified, get it from original model
  if (is.null(var.sb)) var.sb <- model@kern@varHyp

  # if only var needs to be substituted, no need to call fgpm
  if (all(is.null(ls_s.sb), is.null(ls_f.sb))) {
    # duplicate the original model to build the updated one
    modelup <- model

    # recover R and set up the training self-covariance matrix with the substituting variance
    R <- tcrossprod(model@preMats$L)/model@kern@varHyp
    K.tt <- var.sb * R

    # build preMats and replace them in the model
    L <- t(chol(K.tt))
    LInvY <- backsolve(L, model@sOut, upper.tri = FALSE)
    modelup@preMats <- list(L = L, LInvY = LInvY)

    # update the variance slot
    modelup@kern@varHyp <- var.sb

  } else if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # the model is re-made if this is the last one in the sequence of requested tasks
    if (all(!is.null(ls_f.sb), is.null(ls_s.sb))) {
      modelup <- fgpm(sIn = model@sIn, fIn = model@fIn, sOut = model@sOut, kerType = model@kern@kerType,
                      f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims, f_basType = model@f_proj@basType,
                      var.hyp = var.sb, ls_s.hyp = model@kern@s_lsHyps, ls_f.hyp = ls_f.sb,
                      trace = trace, pbars = pbars, control.optim = control.optim)
    } else if(all(!is.null(ls_s.sb), is.null(ls_f.sb))) {
      modelup <- fgpm(sIn = model@sIn, fIn = model@fIn, sOut = model@sOut, kerType = model@kern@kerType,
                      f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims, f_basType = model@f_proj@basType,
                      var.hyp = var.sb, ls_s.hyp = ls_s.sb, ls_f.hyp = model@kern@f_lsHyps,
                      trace = trace, pbars = pbars, control.optim = control.optim)
    } else {
      modelup <- fgpm(sIn = model@sIn, fIn = model@fIn, sOut = model@sOut, kerType = model@kern@kerType,
                      f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims, f_basType = model@f_proj@basType,
                      var.hyp = var.sb, ls_s.hyp = ls_s.sb, ls_f.hyp = ls_f.sb,
                      trace = trace, pbars = pbars, control.optim = control.optim)
    }

  } else if (model@type == "functional") { # functional-input case *******************************************
    # the model is re-made if this is the last one in the sequence of requested tasks
    modelup <- fgpm(fIn = model@fIn, sOut = model@sOut, kerType = model@kern@kerType,
                    f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims,
                    f_basType = model@f_proj@basType, var.hyp = var.sb, ls_f.hyp = ls_f.sb,
                    trace = trace, pbars = pbars, control.optim = control.optim)

  } else { # scalar-input case *******************************************
    # the model is re-made if this is the last one in the sequence of requested tasks
    modelup <- fgpm(sIn = model@sIn, sOut = model@sOut, kerType = model@kern@kerType,
                    var.hyp = var.sb, ls_s.hyp = ls_s.sb,
                    trace = trace, pbars = pbars, control.optim = control.optim)

  }

  return(modelup)
}
# ----------------------------------------------------------------------------------------------------------


# Function to substitute some Hyperparameters
# ----------------------------------------------------------------------------------------------------------
upd_reeHypers <- function(model, var.re, ls_s.re, ls_f.re, trace, pbars, control.optim) {
  # var is always necessary, so if no required to re-estimate, get it from original model
  if (!isTRUE(var.re)) var.up <- model@kern@varHyp else var.up <- NULL

  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    if (!isTRUE(ls_s.re)) ls_s.up <- model@kern@s_lsHyps else ls_s.up <- NULL
    if (!isTRUE(ls_f.re)) ls_f.up <- model@kern@f_lsHyps else ls_f.up <- NULL

    # the model is always re-made
    modelup <- fgpm(sIn = model@sIn, fIn = model@fIn, sOut = model@sOut, kerType = model@kern@kerType,
                    f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims, f_basType = model@f_proj@basType,
                    var.hyp = var.up, ls_s.hyp = ls_s.up, ls_f.hyp = ls_f.up,
                    trace = trace, pbars = pbars, control.optim = control.optim)

  } else if (model@type == "functional") { # functional-input case *******************************************
    if (!isTRUE(ls_f.re)) ls_f.up <- model@kern@f_lsHyps else ls_f.up <- NULL

    # the model is always re-made
    modelup <- fgpm(fIn = model@fIn, sOut = model@sOut, kerType = model@kern@kerType,
                    f_disType = model@kern@f_disType, f_pdims = model@f_proj@pdims, f_basType = model@f_proj@basType,
                    var.hyp = var.up, ls_s.hyp = ls_s.up, ls_f.hyp = ls_f.up,
                    trace = trace, pbars = pbars, control.optim = control.optim)

  } else { # scalar-input case *******************************************
    if (!isTRUE(ls_s.re)) ls_s.up <- model@kern@s_lsHyps else ls_s.up <- NULL

    # the model is always re-made
    modelup <- fgpm(sIn = model@sIn, sOut = model@sOut, kerType = model@kern@kerType,
                    var.hyp = var.up, ls_s.hyp = ls_s.up, ls_f.hyp = ls_f.up,
                    trace = trace, pbars = pbars, control.optim = control.optim)

  }

  return(modelup)
}
# ----------------------------------------------------------------------------------------------------------
