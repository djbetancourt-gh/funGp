# ==========================================================================================================
# Validators
# ==========================================================================================================

# funGp validator
# ----------------------------------------------------------------------------------------------------------
checkVal_funGp <- function(env){
  if (all(!is.null(env$sIn), !is.null(env$fIn))) { # Hybrid-input case *******************************************
    # consistency in number of points
    if (length(unique(c(nrow(env$sIn), as.numeric(sapply(env$fIn, nrow)), length(env$sOut)))) > 1) {
      stop("Inconsistent number of points. Please check that sIn, sOut and each matrix in fIn have all the same number of rows.")
    }

    # consistency in projection dimension
    if (!is.null(env$fpDims)) {
      if (length(env$fpDims) != length(env$fIn)) {
        stop(paste("Inconsistent number of projection dimensions. The functional input list has", length(env$fIn), "elements, but",
                   length(env$fpDims), "projection dimensions were specified."))
      }
    }

  } else if(!is.null(env$fIn)) { # functional-input case ***************************************
    # check validity and consistency of user inputs
    if (length(unique(c(as.numeric(sapply(env$fIn, nrow)), length(env$sOut)))) > 1) {
      stop("Inconsistent number of points. Please check that sOut and each matrix in fIn have all the same number of rows.")
    }

    # consistency in projection dimension
    if (!is.null(env$fpDims)) {
      if (length(env$fpDims) != length(env$fIn)) {
        stop(paste("Inconsistent number of projection dimensions. The functional input list has", length(env$fIn), "elements, but",
                   length(env$fpDims), "projection dimensions were specified."))
      }
    }

  } else if(!is.null(env$sIn)) { # scalar-input case *******************************************
    # check validity and consistency of user inputs
    if (nrow(env$sIn) != length(env$sOut)) {
      stop("Inconsistent number of points. Please check that sIn and sOut have the same number of rows.")
    }
  }
}
# ----------------------------------------------------------------------------------------------------------


# predict and simulate validator
# ----------------------------------------------------------------------------------------------------------
checkVal_pred_and_sim <- function(env){
  # recover the model
  model <- env$model

  # identify if the check is for prediction or simulation and set data sructures and text fields accordingly
  if (all(is.null(env$sIn.sm), is.null(env$fIn.sm))) { # is for prediction
    if(!is.null(env$sIn.pr)) {sIn.ch <- env$sIn.pr; sName <- "sIn.pr"} else sIn.ch <- NULL
    if(!is.null(env$fIn.pr)) {fIn.ch <- env$fIn.pr; fName <- "fIn.pr"} else fIn.ch <- NULL
    taskName <- "prediction"

  } else { # is for simulation
    if(!is.null(env$sIn.sm)) {sIn.ch <- env$sIn.sm; sName <- "sIn.sm"} else sIn.ch <- NULL
    if(!is.null(env$fIn.sm)) {fIn.ch <- env$fIn.sm; fName <- "fIn.sm"} else fIn.ch <- NULL
    taskName <- "simulation"
  }

  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # consistency in data structures
    if (all(is.null(sIn.ch), is.null(fIn.ch))) {
      stop(paste("Invalid input. The model has both, scalar and functional inputs. Please provide valid scalar and\n ",
                 "functional points to proceed with the ", taskName, "."))
    } else if (is.null(fIn.ch)) {
      stop(paste("Inconsistent data structures. The model has both, scalar and functional inputs, but only scalar points\n ",
                 "for ", taskName, " were specified. Please provide also valid new functional points to proceed with the ", taskName, "."))
    } else if (is.null(sIn.ch)) {
      stop(paste("Inconsistent data structures. The model has both, scalar and functional inputs, but only functional points\n ",
                 "for ", taskName, " were specified. Please provide also valid new scalar points to proceed with the ", taskName, "."))
    }

    # consistency in number of points
    if (!all(nrow(sIn.ch) == sapply(fIn.ch, nrow))) {
      stop("Inconsistent number of points. Please check that ", sName, " and each matrix in ", fName, " have all the same number of rows.")
    }

  } else if(model@type == "functional") { # functional-input case ***************************************
    # consistency in data structures
    if (all(!is.null(sIn.ch), !is.null(sIn.ch))) {
      stop(paste("Inconsistent data structures. The model has only functional inputs, but also scalar points\n ",
                 "for ", taskName, " were specified. Please provide only new functional points to proceed with the ", taskName, "."))
    }
    if (!is.null(sIn.ch)) {
      stop(paste("Inconsistent data structures. The model has only functional inputs, but scalar points\n ",
                 "for ", taskName, " were specified. Please provide new functional points instead to proceed with the ", taskName, "."))
    }
    if (is.null(fIn.ch)) {
      stop(paste("Invalid input. The model has functional inputs. Please provide valid new functional points to proceed with\n ",
                 "the ", taskName, "."))
    }

  } else { # scalar-input case *******************************************
    # consistency in data structures
    if (all(!is.null(sIn.ch), !is.null(fIn.ch))) {
      stop(paste("Inconsistent data structures. The model has only scalar inputs, but also functional points\n ",
                 "for ", taskName, " were specified. Please provide only new scalar points to proceed with the ", taskName, "."))
    }
    if (!is.null(fIn.ch)) {
      stop(paste("Inconsistent data structures. The model has only scalar inputs, but functional points\n ",
                 "for ", taskName, " were specified. Please provide new scalar points instead to proceed with the ", taskName, "."))
    }
    if (is.null(sIn.ch)) {
      stop(paste("Invalid input. The model has scalar inputs. Please provide valid new scalar points to proceed with\n ",
                 "the ", taskName, "."))
    }
  }
}
# ----------------------------------------------------------------------------------------------------------


# upd_del validator
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


# upd_subData validator
# ----------------------------------------------------------------------------------------------------------
check_subData <- function(env) {
  # recover the model
  model <- env$model
  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # consistency in number of points
    if (!check_nrows_mlm(env$sIn.sb, env$fIn.sb, env$sOut.sb)) {
      stop("Inconsistent number of points. Please check that sIn.sb, each matrix in fIn.sb and sOut.sb have all the same number\n of rows.")
    }
    if (nrow(env$sIn.sb) != length(env$ind.sb)) {
      stop(paste("Inconsistent number of points in your replacement index vector ind.sb. Please check that it has the same number of\n ",
                 "rows than sIn.sb, each matrix in fIn.sb and sOut.sb.", sep = ""))
    }

  } else if (model@type == "functional") { # functional-input case *******************************************
    # consistency in number of points
    if (!check_nrows_lm(env$fIn.sb, env$sOut.sb)) {
      stop("Inconsistent number of points. Please check that each matrix in fIn.sb and sOut.sb have all the same number of rows.")
    }
    if (nrow(env$fIn.sb[[1]]) != length(env$ind.sb)) {
      stop(paste("Inconsistent number of points in your replacement index vector ind.sb. Please check that it has the same number of\n ",
                 "rows than each matrix in fIn.sb and sOut.sb.", sep = ""))
    }

  } else { # scalar-input case *******************************************
    # consistency in number of points
    if (!check_nrows_mm(env$sIn.sb, env$sOut.sb)) {
      stop("Inconsistent number of points. Please check that sIn.sb and sOut.sb have the same number of rows.")
    }
    if (nrow(env$sIn.sb) != length(env$ind.sb)) {
      stop(paste("Inconsistent number of points in your replacement index vector ind.sb. Please check that it has the same number of\n ",
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


# upd_add validator
# ----------------------------------------------------------------------------------------------------------
check_add <- function(env) {
  # recover the model
  model <- env$model

  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # consistency in number of points
    if (!check_nrows_mlm(env$sIn.nw, env$fIn.nw, env$sOut.nw)) {
      stop("Inconsistent number of points. Please check that sIn.nw, each matrix in fIn.nw and sOut.nw have all the same number\n of rows.")
    }

  } else if (model@type == "functional") { # functional-input case *******************************************
    # consistency in number of points
    if (!check_nrows_lm(env$fIn.nw, env$sOut.nw)) {
      stop("Inconsistent number of points. Please check that each matrix in fIn.nw and sOut.nw have all the same number of rows.")
    }

  } else { # scalar-input case *******************************************
    # consistency in number of points
    if (!check_nrows_mm(env$sIn.nw, env$sOut.nw)) {
      stop("Inconsistent number of points. Please check that sIn.nw and sOut.nw have the same number of rows.")
    }
  }
}
# ----------------------------------------------------------------------------------------------------------

# upd_subHypers validator
# ----------------------------------------------------------------------------------------------------------
check_subHypers <- function(env) {
  # recover the model
  model <- env$model

  # validity of variance if provided: should be > 0
  if (!is.null(env$var.sb)) if (env$var.sb <= 0) stop("The variance should be a positive real number. Please check your inputs.")

  if (model@type == "hybrid") { # Hybrid-input case *******************************************
    # consistency of dimension of scalar length-scale vector if provided
    if (!is.null(env$ls_s.sb)) {
      if (length(model@kern@s_lsHyps) != length(env$ls_s.sb)) {
        stop(paste("The model has ", length(model@kern@s_lsHyps), " scalar length-scale parameters, but you provided",
                   length(env$ls_s.sb), " instead for substitution. Please check your inputs.", sep = ""))
      }
      if (env$ls_s.sb <= 0) stop("Length-scale parameters should be positive real numbers. Please check your ls_s.sb vector.")
    }

    # consistency of dimension of functional length-scale vector if provided
    if (!is.null(env$ls_f.sb)) {
      if (length(model@kern@f_lsHyps) != length(env$ls_f.sb)) {
        stop(paste("The model has ", length(model@kern@f_lsHyps), " functional length-scale parameters, but you provided",
                   length(env$ls_f.sb), " instead for substitution. Please check your inputs.", sep = ""))
      }
      if (env$ls_f.sb <= 0) stop("Length-scale parameters should be positive real numbers. Please check your ls_f.sb vector.")
    }


  } else if (model@type == "functional") { # functional-input case *******************************************
    # consistency of dimension of functional length-scale vector if provided
    if (!is.null(env$ls_f.sb)) {
      if (length(model@kern@f_lsHyps) != length(env$ls_f.sb)) {
        stop(paste("The model has ", length(model@kern@f_lsHyps), " functional length-scale parameters, but you provided",
                   length(env$ls_f.sb), " instead for substitution. Please check your inputs.", sep = ""))
      }
      if (env$ls_f.sb <= 0) stop("Length-scale parameters should be positive real numbers. Please check your ls_f.sb vector.")
    }

  } else { # scalar-input case *******************************************
    # consistency of dimension of scalar length-scale vector if provided
    if (!is.null(env$ls_s.sb)) {
      if (length(model@kern@s_lsHyps) != length(env$ls_s.sb)) {
        stop(paste("The model has ", length(model@kern@s_lsHyps), " scalar length-scale parameters, but you provided",
                   length(env$ls_s.sb), " instead for substitution. Please check your inputs.", sep = ""))
      }
      if (env$ls_s.sb <= 0) stop("Length-scale parameters should be positive real numbers. Please check your ls_s.sb vector.")
    }
  }
}
# ----------------------------------------------------------------------------------------------------------



# ==========================================================================================================
# Checkers
# ==========================================================================================================

# checks if the number of rows in two matrices M1 and M2, and one list L match
# ----------------------------------------------------------------------------------------------------------
check_nrows_mlm <- function(M1, L, M2){
  return(all(nrow(M1) == c(sapply(L, nrow), nrow(M2))))
}
# ----------------------------------------------------------------------------------------------------------


# checks if the number of rows in one list L and one matrix M match
# ----------------------------------------------------------------------------------------------------------
check_nrows_lm <- function(L, M){
  return(all(nrow(M) == sapply(L, nrow)))
}
# ----------------------------------------------------------------------------------------------------------


# checks if the number of rows in two matrices M1 and M2 match
# ----------------------------------------------------------------------------------------------------------
check_nrows_mm <- function(M1, M2){
  return(nrow(M1) == nrow(M2))
}
# ----------------------------------------------------------------------------------------------------------


# checks if there is any duplicate when comparing the duple (sCand, fCand) with the duple (sBench, fBench)
# ----------------------------------------------------------------------------------------------------------
check_duplicates_SF <- function(sBench, fBench, sCand, fCand){
  # merge the benchmark inputs into a single matrix
  bchM <- cbind(sBench, do.call(cbind, fBench))

  # merge the candidate inputs into a single matrix
  cndM <- cbind(sCand, do.call(cbind, fCand))

  # identify duplicates
  if (isTRUE(all.equal(bchM, cndM))) {
    ind.dp <- which(duplicated(bchM) | duplicated(bchM[nrow(bchM):1, ])[nrow(bchM):1])
  } else {
    ind.dp <- which(tail(duplicated(rbind(bchM, cndM)), nrow(cndM)))
  }

  return(ind.dp)
}
# ----------------------------------------------------------------------------------------------------------


# checks if there is any duplicate when comparing fCand with fBench
# ----------------------------------------------------------------------------------------------------------
check_duplicates_F <- function(fBench, fCand){
  # merge the benchmark inputs into a single matrix
  bchM <- do.call(cbind, fBench)

  # merge the candidate inputs into a single matrix
  cndM <- do.call(cbind, fCand)

  # identify duplicates
  if (isTRUE(all.equal(bchM, cndM))) {
    ind.dp <- which(duplicated(bchM) | duplicated(bchM[nrow(bchM):1, ])[nrow(bchM):1])
  } else {
    ind.dp <- which(tail(duplicated(rbind(bchM, cndM)), nrow(cndM)))
  }

  return(ind.dp)
}
# ----------------------------------------------------------------------------------------------------------


# checks if there is any duplicate when comparing sCand with sBench
# ----------------------------------------------------------------------------------------------------------
check_duplicates_S <- function(sBench, sCand, oCand, iCand){
  # merge the benchmark inputs into a single matrix
  bchM <- sBench

  # merge the candidate inputs into a single matrix
  cndM <- sCand

  # identify duplicates
  if (isTRUE(all.equal(bchM, cndM))) {
    ind.dp <- which(duplicated(bchM) | duplicated(bchM[nrow(bchM):1, ])[nrow(bchM):1])
  } else {
    ind.dp <- which(tail(duplicated(rbind(bchM, cndM)), nrow(cndM)))
  }

  return(ind.dp)
}
# ----------------------------------------------------------------------------------------------------------
