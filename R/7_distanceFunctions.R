# Dispatcher for computation of scalar matrices
# -------------------------------------------------------------------------------------------------------------------------------------
setDistMatrix_F <- function(fIn1, fIn2, J, disType){ # fIn1, fIn2 and J are lists. disType is a vector
  Dl <- list()
  for (i in 1:length(J)) {
    switch(disType[i],
           "L2_byindex" = {
             Dl <- c(Dl, L2_byindex(fIn1[[i]], fIn2[[i]]))
           },

           "L2_bygroup" = {
             # Dl[[i]] <- L2_bygroup(fIn1[[i]], fIn2[[i]], J[[i]])
             Dl <- c(Dl, list(L2_bygroup(fIn1[[i]], fIn2[[i]], J[[i]])))
           })
  }
  return(Dl)
}
# -------------------------------------------------------------------------------------------------------------------------------------

setDistMatrix_S <- function(sIn1, sIn2){ # fIn1, fIn2 and J are lists.
  return(L2_byindex(sIn1, sIn2))
}


# computes the matrix of distances using the |.| norm for each difference (s1[i,v] - s2[j,v]) for i,j=1,...n, and
# v fixed from {1,...,ds} such that |d| = abs(d), with d = s1[i,v] - s2[j,v]. Performs the same computation for each
# column in s1 and s2, which points out to a given input variable in the model
L2_byindex <- function(s1, s2){
  Dl <- list()
  for (i in 1:ncol(s1)) {
    Dl[[i]] <- abs(outer(s1[,i], s2[,i], '-'))
  }
  return(Dl)
}

# computes the matrix of distances using the |.|J norm for each difference (f1[i,] - f2[j,]) for i,j=1,...n, such that
# |d|J = sqrt(d' J d), with d = f1[i,] - f2[j,]
L2_bygroup <- function(f1, f2, J){
  f1J <- f1 %*% J
  f2J <- f2 %*% J
  Dm <- matrix(0L, nrow = nrow(f1), ncol = nrow(f2))
  for (a in 1:nrow(J)) {
    Dm <- Dm + (outer(f1[,a], f2[,a], '-') * outer(f1J[,a], f2J[,a], '-'))
  }
  return(sqrt(pmax(Dm,0))) # pmax to prevent numerical problems like slightly negative distances
}

getOwners <- function(d, disType, p){
  owners <- c()
  for (i in 1:d) {
    if (disType[i] == "L2_bygroup") {
      owners <- c(owners, i)
    } else {
      owners <- c(owners, rep(i, p[i]))
    }
  }
  return(owners)
}
