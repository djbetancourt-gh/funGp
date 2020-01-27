check_nrows_mlm <- function(M1, L, M2){
  return(all(nrow(M1) == c(sapply(L, nrow), nrow(M2))))
}
# ----------------------------------------------------------------------------------------------------------

check_nrows_lm <- function(L, M){
  return(all(nrow(M) == sapply(L, nrow)))
}
# ----------------------------------------------------------------------------------------------------------

check_nrows_mm <- function(M1, M2){
  return(nrow(M1) == nrow(M2))
}
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

check_duplicates_F <- function(fBench, fCand, oCand, iCand){
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
