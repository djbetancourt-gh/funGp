# ==============================================================================
# Validators
# ==============================================================================

# funGp validator
# ------------------------------------------------------------------------------
checkVal_fgpm <- function(env){
  if (all(!is.null(env$sIn), !is.null(env$fIn))) { # Hybrid-input case *********
    # consistency in number of points
    if (length(unique(c(nrow(env$sIn), as.numeric(sapply(env$fIn, nrow)),
                        length(env$sOut)))) > 1) {
      stop("Inconsistent number of points. Please check that 'sIn', 'sOut' ",
           " and each matrix in 'fIn' have all the same number of rows.")
    }

    # consistency in projection dimension
    if (!is.null(env$f_pdims)) {
      if (length(env$f_pdims) != length(env$fIn)) {
        stop("Inconsistent number of projection dimensions. The ",
             "functional input list has ", length(env$fIn),
             "elements, but ", length(env$f_pdims),
             " projection dimensions were specified.")
      }
    }

    # consistency between length of 'ls_s.hyp' and the number of scalar inputs
    if (!is.null(env$ls_s.hyp)) {
      if (length(env$ls_s.hyp) != ncol(env$sIn)) {
        stop("Inconsistent number of length-scale parameters. ",
             "Please check that the vector 'ls_s.hyp' has as many ",
             "elements as the number of columns in 'sIn', or leave",
             "this argument empty and the scalar length-scale ",
             "parameters will be estimated from the data.")
      }
    }

    # consistency between length of 'ls_f.hyp' and the sum of effective dimensions
    if (!is.null(env$ls_f.hyp)) {
      od <- sapply(env$fIn, ncol) # original dimensions
      pd <- env$f_pdims # coded projection dimensions (including zeros)
      zrs <- which(pd == 0)
      pd[zrs] <- od[zrs] # effective projection dimensions (zeros replaced by original dimension)
      sumdims <- sum(sapply(seq_along(pd),
                            function(i) if(env$f_disType[i] == "L2_bygroup") return(1) else return(pd[i])))
      if (length(env$ls_f.hyp) != sumdims) {
        stop("Inconsistent number of length-scale parameters. ",
             "Please check that the vector 'ls_f.hyp' has as many ",
             "elements as effective dimensions.\n",
             "Consider checking ?fgpm for details. You can also leave",
             "this argument empty and the functional length-scale",
             "parameters will be  estimated from the data.")
      }
    }

  } else if (!is.null(env$fIn)) { # functional-input case **********************
    # check validity and consistency of user inputs
    if (length(unique(c(as.numeric(sapply(env$fIn, nrow)), length(env$sOut)))) > 1) {
      stop("Inconsistent number of points. Please check that 'sOut' ",
           "and each matrix in 'fIn' have all the same number of rows.")
    }

    # consistency in projection dimension
    if (!is.null(env$f_pdims)) {
      if (length(env$f_pdims) != length(env$fIn)) {
        if (length(env$f_pdims) == 1) {
          stop("Inconsistent number of projection dimensions. ",
               "The functional input list has ", length(env$fIn), " elements, ",
               "but only ", length(env$f_pdims), " projection dimension was ",
               "specified.")
        } else {
          stop("Inconsistent number of projection dimensions. The functional ",
               "input list has ", length(env$fIn), " elements, but ",
               length(env$f_pdims), " projection dimensions were ",
               "specified.")
        }
      }
    }

    # consistency between length of ls_f.hyp and the sum of effective dimensions
    if (!is.null(env$ls_f.hyp)) {
      od <- sapply(env$fIn, ncol) # original dimensions
      pd <- env$f_pdims # coded projection dimensions (including zeros)
      zrs <- which(pd == 0)
      pd[zrs] <- od[zrs] # effective projection dimensions (zeros replaced by original dimension)
      sumdims <- sum(sapply(seq_along(pd),
                            function(i) if(env$f_disType[i] == "L2_bygroup") return(1) else return(pd[i])))
      if (length(env$ls_f.hyp) != sumdims) {
        stop("Inconsistent number of length-scale parameters. Please check that ",
             "the vector 'ls_f.hyp' has as many elements as effective ",
             "dimensions.\n  Consider checking '?fgpm' for details. You can ",
             "also leave this argument empty and the functional length-scale",
             "parameters will be\n  estimated from the data.")
      }
    }

  } else if(!is.null(env$sIn)) { # scalar-input case ************************
    # check validity and consistency of user inputs
    if (nrow(env$sIn) != length(env$sOut)) {
      stop("Inconsistent number of points. Please check that 'sIn' ",
           "and 'sOut' have the same number of rows.")
    }

    if (!is.null(env$ls_s.hyp)) {
      if (length(env$ls_s.hyp) != ncol(env$sIn)) {
        stop("Inconsistent number of length-scale parameters. Please check ",
             "that the vector 'ls_s.hyp' has as many elements as the ",
             "number\n of columns in 'sIn', or leave this argument empty and ",
             "the scalar length-scale parameters will be estimated from the ",
             "data.")
      }
    }
  }

  # validity of the nugget
  if (!is.numeric(env$nugget)) {
    stop("The nugget should be numeric.")
  } else if (env$nugget < 0) {
    stop("The nugget should be a nonnegative number.")
  }

  # validity of n.starts
  if (!is.numeric(env$n.starts)) {
    stop("The argument n.starts should be numeric.")
  } else if (!check.int(env$n.starts)) {
    stop("The argument n.starts should be an integer value.")
  } else if (env$n.starts <= 0) {
    stop("The argument n.starts should be a positive integer.")
  }

  # validity of n.presample
  if (!is.numeric(env$n.presample)) {
    stop("The argument n.presample should be numeric.")
  } else if (!check.int(env$n.presample)) {
    stop("The argument n.presample should be an integer value.")
  } else if (env$n.presample <= 0) {
    stop("The argument n.presample should be a positive integer.")
  }
}
# -----------------------------------------------------------------------------


# predict and simulate validator
# -----------------------------------------------------------------------------
checkVal_pred_and_sim <- function(env) {
  # recover the model
  model <- env$model

  # identify if the check is for prediction or simulation and set
  # data sructures and text fields accordingly
  if (all(is.null(env$sIn.sm), is.null(env$fIn.sm))) { # is for prediction
    if (!is.null(env$sIn.pr)) {
      sIn.ch <- env$sIn.pr
      sName <- "sIn.pr"
    } else sIn.ch <- NULL
    if (!is.null(env$fIn.pr)) {
      fIn.ch <- env$fIn.pr
      fName <- "fIn.pr"
    } else fIn.ch <- NULL
    taskName <- "prediction"

  } else { # is for simulation
    if (!is.null(env$sIn.sm)) {
      sIn.ch <- env$sIn.sm
      sName <- "sIn.sm"
    } else sIn.ch <- NULL
    if (!is.null(env$fIn.sm)) {
      fIn.ch <- env$fIn.sm
      fName <- "fIn.sm"
    } else fIn.ch <- NULL
    taskName <- "simulation"
  }

  if (model@type == "hybrid") { # Hybrid-input case ****************************
    # consistency in data structures
    if (all(is.null(sIn.ch), is.null(fIn.ch))) {
      stop("Invalid input. The model has both, scalar and functional inputs. ",
           "Please provide valid scalar and\n functional points to proceed ",
           "with the ", taskName, ".")
    } else if (is.null(fIn.ch)) {
      stop("Inconsistent data structures. The model has both, scalar and ",
           "functional inputs, but only scalar points\n for ",
           taskName, " were specified. Please provide also valid new ",
           "functional points to proceed with the ", taskName, ".")
    } else if (is.null(sIn.ch)) {
      stop("Inconsistent data structures. The model has both, scalar and ",
           "functional inputs, but only functional points\n for ",
           taskName, " were specified. Please provide also valid new scalar",
           " points to proceed with the ", taskName, ".")
    }

    # consistency in number of points
    if (!all(nrow(sIn.ch) == sapply(fIn.ch, nrow))) {
      stop("Inconsistent number of points. Please check that ", sName, " and ",
           "each matrix in ", fName, " have all the same number of rows.")
    }

  } else if(model@type == "functional") { # functional-input case **************
    # consistency in data structures
    if (all(!is.null(sIn.ch), !is.null(sIn.ch))) {
      stop("Inconsistent data structures. The model has only functional ",
           "inputs, but also scalar points\n for ", taskName, " were ",
           "specified. Please provide only new functional points to ",
           "proceed with the ", taskName, ".")
    }
    if (!is.null(sIn.ch)) {
      stop("Inconsistent data structures. The model has only functional ",
           "inputs, but scalar points\n for ", taskName, " were ",
           "specified. Please provide new functional points instead ",
           "to proceed with the ", taskName, ".")
    }
    if (is.null(fIn.ch)) {
      stop("Invalid input. The model has functional inputs. Please ",
           "provide valid new functional points to proceed with\n the ",
           taskName, ".")
    }

  } else { # scalar-input case *************************************************
    # consistency in data structures
    if (all(!is.null(sIn.ch), !is.null(fIn.ch))) {
      stop("Inconsistent data structures. The model has only scalar inputs, ",
           "but also functional points\n for ", taskName, " were specified. ",
           "Please provide only new scalar points to proceed with the ",
           taskName, ".")
    }
    if (!is.null(fIn.ch)) {
      stop("Inconsistent data structures. The model has only scalar inputs, ",
           "but functional points\n for ", taskName, " were specified. ",
           "Please provide new scalar points instead to proceed with the ",
           taskName, ".")
    }
    if (is.null(sIn.ch)) {
      stop("Invalid input. The model has scalar inputs. Please provide valid ",
           "new scalar points to proceed with\n the ", taskName, ".")
    }
  }
}
# ------------------------------------------------------------------------------


# upd_del validator
# ------------------------------------------------------------------------------
check_del <- function(env) {
  # recover the model
  model <- env$model

  # is there any duplicate?
  if (any(duplicated(env$ind.dl))) {
    warning("Some elements in the deletion index vector 'ind.dl' are ",
            "duplicated. Duplicates are dropped.")
    env$ind.dl <- unique(env$ind.dl)
  }

  # are indices in the correct range?
  if (any(env$ind.dl < 1) || any(env$ind.dl > model@n.tot)) {
    stop("Some elements in the detelion index vector 'ind.dl' are not in ",
         "'[1, model@n.tot]'. Please check your vector.")
  }

  # is the maximum number of points to delete respected?
  if (model@n.tot - length(env$ind.dl) < 2) {
    stop("The number of points to delete cannot exceed 'model@n.tot - 2'. ",
         "Please check your vector.")
  }
  return(env$ind.dl)
}
# ------------------------------------------------------------------------------


# upd_subData validator
# ------------------------------------------------------------------------------
check_subData <- function(env) {
  # recover the model
  model <- env$model
  if (model@type == "hybrid") { # Hybrid-input case ****************************
    # consistency in number of points
    if (!check_nrows_mlm(env$sIn.sb, env$fIn.sb, env$sOut.sb)) {
      stop("Inconsistent number of points. Please check that 'sIn.sb', each ",
           "matrix in 'fIn.sb' and 'sOut.sb' have all the same number\n'",
           "of rows.")
    }
    if (nrow(env$sIn.sb) != length(env$ind.sb)) {
      stop("Inconsistent number of points in your replacement index vector ",
           "'ind.sb'. Please check that it has the same number of\n ",
           "rows as 'sIn.sb', each matrix in 'fIn.sb' and 'sOut.sb'.")
    }

  } else if (model@type == "functional") { # functional-input case *************
    # consistency in number of points
    if (!check_nrows_lm(env$fIn.sb, env$sOut.sb)) {
      stop("Inconsistent number of points. Please check that all matrices ",
           "in 'fIn.sb' and 'sOut.sb' have all the same number of rows.")
    }
    if (nrow(env$fIn.sb[[1]]) != length(env$ind.sb)) {
      stop("Inconsistent number of points in your replacement index vector ",
           "'ind.sb'. Please check that it has the same number of\n ",
           "rows than each matrix in fIn.sb and sOut.sb.")
    }

  } else { # scalar-input case *************************************************
    # consistency in number of points
    if (!check_nrows_mm(env$sIn.sb, env$sOut.sb)) {
      stop("Inconsistent number of points. Please check that 'sIn.sb' and ",
           "'sOut.sb' have the same number of rows.")
    }
    if (nrow(env$sIn.sb) != length(env$ind.sb)) {
      stop("Inconsistent number of points in your replacement index vector ",
           "'ind.sb'. Please check that it has the same number of\n ",
           "rows than sIn.sb and sOut.sb.")
    }
  }

  # validity of subtituting index
  # should be an integer
  if (any(env$ind.sb %% 1 != 0)) {
    stop("Replacement indices should be integer numbers. Please check the ",
         "following positions of your 'ind.sb' vector: \n",
         paste(which(env$ind.sb %% 1 != 0), collapse = ", "),
         ".")
  }
  # should be in [1, n.tot]
  if (any(env$ind.dl < 1) || any(env$ind.dl > model@n.tot)) {
    stop("Replacement indices should be integers in '[1, model@n.tot]'. ",
         "Please check the following positions of your 'ind.sb' vector: \n",
         paste(which(!(env$ind.sb > 0 & env$ind.sb < model@n.tr)),
               collapse = ", "),
         ".")
  }
  # should not contain duplicates
  if (anyDuplicated(env$ind.sb)) {
    stop("Replacement indices should be unique. Please check the following ",
         "positions of your ind.sb vector: ",
         paste(which(duplicated(env$ind.sb) | duplicated(env$ind.sb, fromLast = TRUE)),
               collapse = " "), ".")
  }
}
# ------------------------------------------------------------------------------


# upd_add validator
# ------------------------------------------------------------------------------
check_add <- function(env) {
  # recover the model
  model <- env$model

  if (model@type == "hybrid") { # Hybrid-input case ****************************
    # consistency in number of points
    if (!check_nrows_mlm(env$sIn.nw, env$fIn.nw, env$sOut.nw)) {
      stop("Inconsistent number of points. Please check that 'sIn.nw', ",
           "each matrix in fIn.nw and sOut.nw have all the same number\n ",
           "of rows.")
    }

  } else if (model@type == "functional") { # functional-input case *************
    # consistency in number of points
    if (!check_nrows_lm(env$fIn.nw, env$sOut.nw)) {
      stop("Inconsistent number of points. Please check that each matrix ",
           "in 'fIn.nw' and 'sOut.nw' have all the same number of rows.")
    }

  } else { # scalar-input case *************************************************
    # consistency in number of points
    if (!check_nrows_mm(env$sIn.nw, env$sOut.nw)) {
      stop("Inconsistent number of points. Please check that 'sIn.nw' ",
           "and 'sOut.nw' have the same number of rows.")
    }
  }
}
# ------------------------------------------------------------------------------


# upd_subHypers validator
# ------------------------------------------------------------------------------
check_subHypers <- function(env) {
  # recover the model
  model <- env$model

  # validity of variance if provided: should be > 0
  if (!is.null(env$var.sb))
    if (env$var.sb <= 0)
      stop("The variance should be a positive real number. ",
           "Please check your inputs.")

  if (model@type == "hybrid") { # Hybrid-input case ****************************
    # consistency of dimension of scalar length-scale vector if provided
    if (!is.null(env$ls_s.sb)) {
      if (length(model@kern@s_lsHyps) != length(env$ls_s.sb)) {
        stop("The model has ", length(model@kern@s_lsHyps), " scalar ",
             "length-scale parameters, but you provided",
             length(env$ls_s.sb), " instead for substitution. Please ",
             "check your inputs.")
      }
      if (any(env$ls_s.sb <= 0))
        stop("Length-scale parameters should be positive real numbers. ",
             "Please check your ls_s.sb vector.")
    }

    # consistency of dimension of functional length-scale vector if provided
    if (!is.null(env$ls_f.sb)) {
      if (length(model@kern@f_lsHyps) != length(env$ls_f.sb)) {
        stop("The model has ", length(model@kern@f_lsHyps), " functional ",
             "length-scale parameters, but you provided ",
             length(env$ls_f.sb), " instead for substitution. Please ",
             "check your inputs.")
      }
      if (any(env$ls_f.sb <= 0))
        stop("Length-scale parameters should be positive real numbers. ",
             "Please check your 'ls_f.sb' vector.")
    }


  } else if (model@type == "functional") { # functional-input case *************
    # consistency of dimension of functional length-scale vector if provided
    if (!is.null(env$ls_f.sb)) {
      if (length(model@kern@f_lsHyps) != length(env$ls_f.sb)) {
        stop("The model has ", length(model@kern@f_lsHyps), " functional ",
             "length-scale parameters, but you provided ",
             length(env$ls_f.sb), " instead for substitution. Please ",
             "check your inputs.")
      }
      if (env$ls_f.sb <= 0)
        stop("Length-scale parameters should be positive real numbers. ",
             "Please check your ls_f.sb vector.")
    }

  } else { # scalar-input case *************************************************
    # consistency of dimension of scalar length-scale vector if provided
    if (!is.null(env$ls_s.sb)) {
      if (length(model@kern@s_lsHyps) != length(env$ls_s.sb)) {
        stop("The model has ", length(model@kern@s_lsHyps), " scalar ",
             "length-scale parameters, but you provided ",
             length(env$ls_s.sb), " instead for substitution. Please ",
             "check your inputs.")
      }
      if (env$ls_s.sb <= 0)
        stop("Length-scale parameters should be positive real numbers. ",
             "Please check your 'ls_s.sb' vector.")
    }
  }
}
# ------------------------------------------------------------------------------



# =============================================================================
# Checkers
# =============================================================================

# checks if the number of rows in two matrices M1 and M2, and one list L match
# ------------------------------------------------------------------------------
check_nrows_mlm <- function(M1, L, M2){
  return(all(nrow(M1) == c(sapply(L, nrow), nrow(M2))))
}
# ------------------------------------------------------------------------------


# checks if the number of rows in one list L and one matrix M match
# ------------------------------------------------------------------------------
check_nrows_lm <- function(L, M){
  return(all(nrow(M) == sapply(L, nrow)))
}
# ------------------------------------------------------------------------------


# checks if the number of rows in two matrices M1 and M2 match
# ------------------------------------------------------------------------------
check_nrows_mm <- function(M1, M2){
  return(nrow(M1) == nrow(M2))
}
# ------------------------------------------------------------------------------


# checks if there is any duplicate when comparing the duple (sCand, fCand) with
# the duple (sBench, fBench)
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------


# checks if there is any duplicate when comparing 'fCand' with 'fBench'
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------


# checks if there is any duplicate when comparing 'sCand' with 'sBench'
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
#
#
#  check_new_inputs
# ------------------------------------------------------------------------------
##' Check the consistency of provided data describing some scalar and
##' functional inputs, possibly in relation with a \code{fgpm}
##' object. This function is to be used when predicting and simulating
##' from a \code{fgpm} object or to create such an object.
##'
##' Note that when passed to \code{newsIn} or \code{newfIn}, any
##' object with length zero will be accepted to specify
##' "no inputs". So a logical or character vector of length zero can
##' be accepted, resulting either in a matrix with zero columns
##' (scalar inputs) or a list with length zero (functional
##' inputs). Also note that numeric vectors are not accepted as
##' one-column matrices for the functional inputs because a functional
##' input should have at least one "time". So some care is needed when
##' a \code{fgpm} object is programmatically created if a dimension
##' reduction is used possibly ending in a one-dimensional central
##' space: the dimension should not be dropped.
##'
##' @title Check the consistency of inputs
##'
##' @param object Optional \code{fgpm} object. If missing, the
##'     function will only check the consistency between \code{newsIn}
##'     and \code{newfIn} as required when creating a new \code{fgpm}
##'     object. If \code{object} is given it must be of class
##'     \code{fgpm} and the consistency of \code{newsIn} and
##'     \code{newfIn} with this object is further checked.
##' @param newsIn A numeric matrix of scalar inputs or an object with
##'     length zero (such as \code{NULL}) if no scalar inputs are to
##'     be used.
##' @param newfIn A list describing the functional inputs or an object
##'     with length zero (such as \code{NULL}) if no functional inputs
##'     are to be used. If some functional inputs are used, the
##'     elements of \code{newfIn} must be numeric matrices with the
##'     same number of rows.
##'
##' @return A list with the following elements
##' \itemize{
##'    \item{n }{The number of observations.}
##'    \item{newsIn }{ The checked and possibly corrected matrix of
##'       scalar inputs. This will always be a numeric matrix with \code{n}
##'       rows, possibly with zero columns.}
##'    \item{newfIn }{The checked and possibly corrected list of functional
##'       inputs. This can be a list with zero elements, but if there are
##'       some elements all will be numeric matrices with the same number
##'       of columns.}
##' }
##' So the "official" way to query about the number of scalar and functional
##' inputs is to use \code{ncol(res$newsIn)} and \code{length(res$newfIn)}
##' where \code{res} is the returned object.
##'
##' @section Caution: There should be at least one input. Should/could
##'     the function also check for names? \bold{This function is
##'     exported only on a temporary basis}, to make the help and
##'     examples visible.
##'
##' @export
##'
##' @examples
##' ## works
##' check_new_inputs(newsIn = runif(4))
##' check_new_inputs(newsIn = runif(4),
##'                  newfIn = list(matrix(runif(24), nrow = 4, ncol = 6)))
##' check_new_inputs(newsIn = matrix(runif(8), nrow = 4, ncol = 2),
##'                  newfIn = list(matrix(runif(24), nrow = 4, ncol = 6)))
##' check_new_inputs(newfIn = list(matrix(runif(24), nrow = 4, ncol = 6)))
##'
##' ## errors in the functional part
##' try(check_new_inputs(newsIn = runif(4), newfIn = runif(4)))
##' try(check_new_inputs(newsIn = runif(4), newfIn = list(runif(4))))
##' try(check_new_inputs(newsIn = runif(4),
##'                      newfIn = list(matrix(runif(5), nrow = 5, ncol = 1))))
##' try(check_new_inputs(newsIn = letters[1:5],
##'                      newfIn = list(matrix(runif(5), nrow = 5, ncol = 1))))
##' try(check_new_inputs(newsIn = runif(5),
##'                      newfIn = list(matrix(letters[1:5], nrow = 5, ncol = 1))))
##' ## Make some test 'fgpm' objects as in 'example("fgpm")'
##' set.seed(100); n.tr <- 25
##' sIn <- expand.grid(x1 = seq(0, 1, length = sqrt(n.tr)),
##'                    x2 = seq(0, 1, length = sqrt(n.tr)))
##' fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
##'             f2 = matrix(runif(n.tr * 22), ncol = 22))
##' sOut <- fgp_BB3(sIn, fIn, n.tr)
##' ms <- fgpm(sIn = sIn, sOut = sOut)
##' mf <- fgpm(fIn = fIn, sOut = sOut)
##' msf <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
##' ## check with the inputs used at creation.
##' res_s <- check_new_inputs(object = ms, newsIn = sIn)
##' res_f <- check_new_inputs(object = mf, newfIn = fIn)
##' res_sf <- check_new_inputs(object = msf, newsIn = sIn, newfIn = fIn)
##'
check_new_inputs <- function(object, newsIn = NULL, newfIn = NULL) {

    ## We have to guess the number 'n' of observations, which in
    ## general may require inspecting 'newsIn' AND 'newfIn'. So
    ## temporary values 'ns' and 'nf' are used.

    ## Check 'newsIn', allowing for a numeric vector
    if (!length(newsIn)) {
        newsIn <- matrix(numeric(0), ncol = 0)
        nFroms <- NA
    } else {
        if (is.data.frame(newsIn)) {
            ## warning("coerce 'newsIn' from \"data.frame\" into \"matrix\". ",
            ##         "The order of the columns matter!")
        }
        newsIn <- as.matrix(newsIn)
        if (!is.numeric(newsIn)) stop("'newsIn' must be numeric")
        nFroms <- nrow(newsIn)
    }

    ## Check 'newfIn' Begin with the case of an object with length 0
    ## which covers the case when 'newfIn' is NULL
    if (!length(newfIn)) {
        newfIn <- list()
        nFromf <- NA
    } else {
        if (!is.list(newfIn) || !all(sapply(newfIn, is.matrix)) ||
            !all(sapply(newfIn, is.numeric)) ||
            (length(unique(sapply(newfIn, nrow))) != 1))
            stop("'newfIn' must be a list of numeric matrices with the same ",
                 "number of rows")
        nFromf <- nrow(newfIn[[1]])
    }

    ## Check 'newfIn' Begin with the case of an object with length 0
    ## which covers the case when 'newfIn' is NULL
    if (!is.na(nFroms) && !is.na(nFromf)) {
        if (nFroms != nFromf)
            stop("'newsIn' and 'newfIn' mismatch. Please check the ",
                 "numbers of rows.")
        n <- nFroms
    } else {
        n <- min(c(nFroms, nFromf), na.rm = TRUE)
        ## force the number of rows of 'newsIn' (note that `nrow<-`
        ## does not exist)
        if (!ncol(newsIn)) dim(newsIn) <- c(n , 0)
    }

    ## Now check the consistency of 'newsIn' and 'newfIn' with 'object' if
    ## needed.
    if (!missing(object)) {
        if (!inherits(object, "fgpm"))
            stop("'object' must inherit from \"fgpm\"")

        if (object@ds != ncol(newsIn)) {
            stop("'object' requires ", object@ds, " scalar inputs")
        }
        if (object@df != length(newfIn)) {
            stop("'object' requires 'object@df' functional inputs")
        }
    }

    list(n = n, newsIn = newsIn, newfIn = newfIn)

}
