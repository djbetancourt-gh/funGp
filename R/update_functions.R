# ==========================================================================================================
# Model rebuilder
# ==========================================================================================================

# # Function to rebuild a model based on new data and/or parameters
# # ----------------------------------------------------------------------------------------------------------
# refunGp <- function(model, sIn = NULL, fIn = NULL, sOut, n.starts = 1, n.presample = 20){
#   # duplicate the original model to build the updated one
#   modelup <- model
#
#   if (model@type == "hybrid") { # Hybrid-input case *******************************************
#     # extract information from original model
#     ds <- model@ds
#     df <- model@df
#     doProj <- model@proj@doProj
#     fpDims <- model@proj@fpDims
#
#     if (doProj) {
#       if (is.null(fpDims)) {
#         fpDims <- rep(3, df)
#       }
#
#       # project functional inputs
#       basis <- fpIn <- J <- list()
#       for (i in 1:df) {
#         if (fpDims[i] > 0) {
#           B <- (eigen(cov(fIn[[i]]))$vectors)[,1:fpDims[i]]
#           fpIn[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn[[i]]))
#           J[[i]] <- t(B) %*% B
#         } else {
#           J[[i]] <- B <- diag(ncol(fIn[[i]]))
#           fpIn[[i]] <- fIn[[i]]
#         }
#         basis[[i]] <- B
#       }
#     } else {
#       fpDims <- rep(0, df)
#       basis <- J <- lapply(fIn, function(m) diag(ncol(m)))
#       fpIn <- fIn
#     }
#
#     # compute scalar distance matrices
#     sMs <- setScalDistance(sIn, sIn)
#
#     # compute functional distance matrices
#     fMs <- setFunDistance(fpIn, fpIn, J)
#
#     # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
#     modelup@preMats <- preMats_SF(sMs, fMs, sOut, model@kern@varHyp, model@kern@s_lsHyps,
#                                   model@kern@f_lsHyps, model@kern@kerType)
#
#     # fill funGpProj slots specific to the hybrid-input case
#     modelup@proj@basis <- basis
#     modelup@proj@coefs <- fpIn
#
#     # fill funGp slots specific to the hybrid-input case
#     modelup@sIn <- sIn
#     modelup@fIn <- fIn
#
#   } else if (model@type == "functional") { # functional-input case *******************************************
#     # extract information from original model
#     df <- model@df
#     doProj <- model@proj@doProj
#     fpDims <- model@proj@fpDims
#
#     if (doProj) {
#       if (is.null(fpDims)) {
#         fpDims <- rep(3, df)
#       }
#
#       # project functional inputs
#       basis <- fpIn <- J <- list()
#       for (i in 1:df) {
#         if (fpDims[i] > 0) {
#           B <- (eigen(cov(fIn[[i]]))$vectors)[,1:fpDims[i]]
#           fpIn[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn[[i]]))
#           J[[i]] <- t(B) %*% B
#         } else {
#           J[[i]] <- B <- diag(ncol(fIn[[i]]))
#           fpIn[[i]] <- fIn[[i]]
#         }
#         basis[[i]] <- B
#       }
#     } else {
#       fpDims <- rep(0, df)
#       basis <- J <- lapply(fIn, function(m) diag(ncol(m)))
#       fpIn <- fIn
#     }
#
#     # compute functional distance matrices
#     fMs <- setFunDistance(fpIn, fpIn, J)
#
#     # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
#     modelup@preMats <- preMats_F(fMs, sOut, model@kern@varHyp, model@kern@f_lsHyps, model@kern@kerType)
#
#     # fill funGpProj slots specific to the hybrid-input case
#     modelup@proj@basis <- basis
#     modelup@proj@coefs <- fpIn
#
#     # fill funGp slots specific to the functional-input case
#     modelup@fIn <- fIn
#
#   } else { # scalar-input case *******************************************
#     # extract information from original model
#     ds <- model@ds
#
#     # compute scalar distance matrices
#     sMs <- setScalDistance(sIn, sIn)
#
#     # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
#     modelup@preMats <- preMats_S(sMs, sOut, model@kern@varHyp, model@kern@s_lsHyps, model@kern@kerType)
#
#     # fill funGp slots specific to the scalar-input case
#     modelup@sIn <- sIn
#   }
#
#   # fill general funGpModel slots
#   modelup@sOut <- sOut
#   modelup@n.tot <- length(sOut)
#
#   return(modelup)
# }
# # ----------------------------------------------------------------------------------------------------------



# ==========================================================================================================
# Deletion function and validator
# ==========================================================================================================

# # Function to delete some data
# # ----------------------------------------------------------------------------------------------------------
# update_InOut_dl.funGp <- function(model, ind.dl) {
#   # check for validty of substituting data
#   ind.dl <- checkVal_InOut_dl(as.list(environment()))
#
#   # extract generic information from the model
#   sOut <- model@sOut
#
#   if (model@type == "hybrid") {
#     # extract inputs from original model
#     sIn <- model@sIn
#     fIn <- model@fIn
#
#     # remove points according to deletion indices
#     sIn <- sIn[-ind.dl,,drop = F]
#     fIn <- lapply(fIn, function(M) M[-ind.dl,])
#     sOut <- sOut[-ind.dl,,drop = F]
#
#     # request new model to refunGp
#     modelup <- refunGp(model = model, sIn = sIn, fIn = fIn, sOut = sOut)
#
#   } else if (model@type == "functional") {
#     # extract inputs from original model
#     fIn <- model@fIn
#
#     # remove points according to deletion indices
#     fIn <- lapply(fIn, function(M) M[-ind.dl,])
#     sOut <- sOut[-ind.dl,,drop = F]
#
#     # request new model to refunGp
#     modelup <- refunGp(model = model, fIn = fIn, sOut = sOut)
#
#   } else {
#     # extract inputs from original model
#     sIn <- model@sIn
#
#     # remove points according to deletion indices
#     sIn <- sIn[-ind.dl,,drop = F]
#     sOut <- sOut[-ind.dl,,drop = F]
#
#     # request new model to refunGp
#     modelup <- refunGp(model = model, sIn = sIn, sOut = sOut)
#   }
#
#   return(modelup)
# }
# # ----------------------------------------------------------------------------------------------------------





# ==========================================================================================================
# Substitution function and validator
# ==========================================================================================================

# # Function to substitute some data
# # ----------------------------------------------------------------------------------------------------------
# update_InOut_sb.funGp <- function(model, sIn.sb, fIn.sb, sOut.sb, ind.sb) {
#   browser()
#   # extract generic information from the model
#   sOut <- model@sOut
#
#   # provide substituting output if not specified by the user
#   if(is.null(sOut.sb)) sOut.sb <- sOut[ind.sb,,drop = F]
#
#   if (model@type == "hybrid") { # Hybrid-input case *******************************************
#     # extract inputs from original model
#     sIn <- model@sIn
#     fIn <- model@fIn
#
#     # provide substituting inputs if not specified by the user
#     if(is.null(sIn.sb)) sIn.sb <- sIn[ind.sb,,drop = F]
#     if(is.null(fIn.sb)) fIn.sb <- lapply(fIn, function(M) M[ind.sb,,drop = F])
#
#     # check for validty of substituting data
#     checkVal_InOut_sb(as.list(environment()))
#
#     # check for duplicates in the substituting points
#     ind.dp <- checkDuplicates_SF(sIn.sb, fIn.sb, sIn.sb, fIn.sb, sOut.sb, ind.sb)
#
#     if (length(ind.sb) == length(ind.dp)) {
#       warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
#                     "Substitution is skipped.", sep = ""))
#       return(model)
#     } else if (length(ind.dp) > 0) {
#       warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
#                     "Duplicate substitute points: ", res$ind.dp, sep = ""))
#       sIn.sb <- sIn.sb[-ind.dp,,drop = F]
#       fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = F])
#       Sout.sb <- sOut.sb[-ind.dp,,drop = F]
#       ind.sb <- ins.db[-ind.dp]
#     }
#
#     # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
#     sIn.exsb <- sIn[-ind.sb,]
#     fIn.exsb <- lapply(fIn, function(M) M[-ind.sb,,drop = F])
#     ind.dp <- checkDuplicates_SF(sIn.exsb, fIn.exsb, sIn.sb, fIn.sb, sOut.sb, ind.sb)
#
#     if (length(ind.sb) == length(ind.dp)) {
#       warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
#                     "contained in the model. The model is returned in its original state.", sep = ""))
#       return(model)
#     } else if (length(ind.dp) > 0) {
#       warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
#                     "Duplicate substitute points: ", res$ind.dp, sep = ""))
#       sIn.sb <- sIn.sb[-ind.dp,,drop = F]
#       fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = F])
#       Sout.sb <- sOut.sb[-ind.dp,,drop = F]
#       ind.sb <- ins.db[-ind.dp]
#     }
#
#     # recover inputs and outputs after duplicates check
#     sIn[ind.sb,] <- sIn.sb
#     fIn <- mapply(function(M, x) {M[ind.sb,] <- x; return(M)}, fIn, fIn.sb)
#     sOut[ind.sb,] <- sOut.sb
#
#     # request new model to refunGp
#     modelup <- refunGp(model = model, sIn = sIn, fIn = fIn, sOut = sOut)
#
#   } else if (model@type == "functional") { # functional-input case *******************************************
#     # extract inputs from original model
#     fIn <- model@fIn
#
#     # provide substituting inputs if not specified by the user
#     if(is.null(fIn.sb)) fIn.sb <- lapply(fIn, function(M) M[ind.sb,,drop = F])
#
#     # check for validty of substituting data
#     checkVal_InOut_sb(as.list(environment()))
#
#     # check for duplicates in the substituting points
#     ind.dp <- checkDuplicates_F(fIn.sb, fIn.sb, sOut.sb, ind.sb)
#
#     if (length(ind.sb) == length(ind.dp)) {
#       warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
#                     "Substitution is skipped.", sep = ""))
#       return(model)
#     } else if (length(ind.dp) > 0) {
#       warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
#                     "Duplicate substitute points: ", res$ind.dp, sep = ""))
#       fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = F])
#       Sout.sb <- sOut.sb[-ind.dp,,drop = F]
#       ind.sb <- ins.db[-ind.dp]
#     }
#
#     # check for duplicates bewteen substituting inputs and existing inputs at not substituting rows
#     fIn.exsb <- lapply(fIn, function(M) M[-ind.sb,,drop = F])
#     ind.dp <- checkDuplicates_F(fIn.exsb, fIn.sb, sOut.sb, ind.sb)
#
#     if (length(ind.sb) == length(ind.dp)) {
#       warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
#                     "contained in the model. The model is returned in its original state.", sep = ""))
#       return(model)
#     } else if (length(ind.dp) > 0) {
#       warning(paste("There are some duplicates in the substituting inputs. Those have been ignored.\n",
#                     "Duplicate substitute points: ", res$ind.dp, sep = ""))
#       fIn.sb <- lapply(fIn.sb, function(M) M[-ind.dp,,drop = F])
#       Sout.sb <- sOut.sb[-ind.dp,,drop = F]
#       ind.sb <- ins.db[-ind.dp]
#     }
#
#     # recover inputs and outputs after duplicates check
#     fIn <- mapply(function(M, x) {M[ind.sb,] <- x; return(M)}, fIn, fIn.sb)
#     sOut[ind.sb,] <- sOut.sb
#
#     # request new model to refunGp
#     modelup <- refunGp(model = model, fIn = fIn, sOut = sOut)
#
#   } else { # scalar-input case *******************************************
#     # extract inputs from original model
#     sIn <- model@sIn
#
#     # remove points according to deletion indices
#     sIn <- sIn[-ind.dl,,drop = F]
#     sOut <- sOut[-ind.dl,,drop = F]
#
#     # request new model to refunGp
#     modelup <- refunGp(model = model, sIn = sIn, sOut = sOut)
#   }
#
#   return(modelup)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#   # browser()
#   # duplicate the model for clarity
#   #modelup <- model
#
#   # extract generic information from user inputs
#   #sOut <- model@sOut
#
#   # provide substituting output if not specified by the user
#   #if(is.null(sOut.sb)) sOut.sb <- sOut[ind.sb,,drop = F]
#   #sOut.sb <- as.matrix(sOut.sb)
#
#   # check which type of model it is
#   if (all(model@ds > 0, model@df > 0)) { # Hybrid-input case *******************************************
#     # extract information from user inputs specific to the hybrid-input case
#     #sIn <- model@sIn
#     #fIn <- model@fIn
#
#     # provide substituting inputs if not specified by the user
#     # if(is.null(sIn.sb)) sIn.sb <- sIn[ind.sb,,drop = F]
#     # if(is.null(fIn.sb)) fIn.sb <- lapply(fIn, function(M) M[ind.sb,,drop = F])
#
#     # check for validty of substituting data
#     # checkVal_InOut_sb(as.list(environment()))
#
#     # check for duplicates in the substituting points
#     # res <- checkDuplicates(sBench = sIn.sb, fBench = fIn.sb, sCand = sIn.sb, fCand = fIn.sb, oCand = sOut.sb, iCand = ind.sb)
#
#     # update substituting data and warn if required
#     # if (length(res$iClean) == 0) {
#     #   warning(paste("No substituting points left after checking for duplicates in the substituting input points. ",
#     #                 "The model is returned in its original state.", sep = ""))
#     #   return(model)
#     # } else {
#     #   if (length(res$ind.dp) > 0) {
#     #     warning(paste("There are some duplicates in the substituting inputs. Duplicates have been ignored.\n",
#     #                   "Duplicate substitute points: ", res$ind.dp, sep = ""))
#     #   }
#     #   sIn.sb <- res$sClean
#     #   fIn.sb <- res$fClean
#     #   sOut.sb <- res$oClean
#     #   ind.sb <- res$iClean
#     # }
#
#     # check for duplicated bewteen substituting inputs and existing inputs at not substituting rows
#     # sIn.exsb <- sIn[-ind.sb,]
#     # fIn.exsb <- lapply(fIn, function(M) M[-ind.sb,,drop = F])
#     # res <- checkDuplicates(sBench = sIn.exsb, fBench = fIn.exsb, sCand = sIn.sb, fCand = fIn.sb, oCand = sOut.sb, iCand = ind.sb)
#
#     # update substituting data
#     if (length(res$iClean) == 0) {
#       warning(paste("No substituting points left after cross-checking for duplicates against the inputs already. ",
#                     "contained in the model. The model is returned in its original state.", sep = ""))
#       return(model)
#     } else {
#       if (length(res$ind.dp) > 0) {
#         warning(paste("There are some duplicates in the substituting inputs. Duplicates have been ignored.\n",
#                       "Duplicate substitute points: ", res$ind.dp, sep = ""))
#       }
#       sIn.sb <- res$sClean
#       fIn.sb <- res$fClean
#       sOut.sb <- res$oClean
#       ind.sb <- res$iClean
#
#       # recover inputs and outputs after duplicates check
#       sIn[res$iClean,] <- res$sClean
#       fIn <- mapply(function(M, x) {M[res$iClean,] <- x; return(M)}, fIn, res$fClean)
#       sOut[res$iClean,] <- res$oClean
#
#       # extract information from previous model specific to the hybrid-input case
#       ds <- model@ds
#       df <- model@df
#       doProj <- model@proj@doProj
#       fpDims <- model@proj@fpDims
#
#       # Extend to other possible cases!!!!!!!!!!!!!!!!!!
#       if (doProj) {
#         # project functional inputs
#         basis <- fpIn <- J <- list()
#         for (i in 1:df) {
#           if (fpDims[i] > 0) {
#             B <- (eigen(cov(fIn[[i]]))$vectors)[,1:fpDims[i]]
#             fpIn[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn[[i]]))
#             J[[i]] <- t(B) %*% B
#           } else {
#             J[[i]] <- B <- diag(ncol(fIn[[i]]))
#             fpIn[[i]] <- fIn[[i]]
#           }
#           basis[[i]] <- B
#         }
#       } else {
#         basis <- J <- lapply(fIn, function(m) diag(ncol(m)))
#         fpIn <- fIn
#       }
#       # compute scalar distance matrices
#       sMs <- setScalDistance(sIn, sIn)
#
#       # compute functional distance matrices
#       fMs <- setFunDistance(fpIn, fpIn, J)
#
#       # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
#       modelup@preMats <- preMats_SF(sMs, fMs, sOut, model@kern@varHyp, model@kern@s_lsHyps,
#                                     model@kern@f_lsHyps, model@kern@kerType)
#
#       # fill funGpProj slots specific to the hybrid-input case
#       modelup@proj@basis <- basis
#       modelup@proj@coefs <- fpIn
#
#       # fill funGp slots specific to the hybrid-input case
#       modelup@sIn <- sIn
#       modelup@fIn <- fIn
#     }
#
#   } else if (model@df > 0) { # functional-input case *******************************************
#     print("I'm functional!")
#
#
#   } else { # scalar-input case *******************************************
#
#   }
#
#   # fill general funGpModel slots
#   modelup@sOut <- sOut
#
#   return(modelup)
# }
# # ----------------------------------------------------------------------------------------------------------
