#' @title Computation of distance matrix for scalar inputs
#' @description Precomputes the distance matrices that will be used later inside the \emph{scalar} term of the kernel function
#'              during the optimization of the hyperparmeters of the metamodel.
#'
#' @param sIn1 a first matrix of scalar input values.
#' @param sIn2 a second matrix of scalar input values.
#' @return A list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#'
#' @keywords internal
#'
#' @author José Betancourt
#' @export
setScalDistance <- function(sIn1, sIn2){
  sMs <- list()
  for (l in 1:ncol(sIn1)) {
    sMs[[l]] <- outer(sIn1[,l], sIn2[,l], '-')
  }
  return(sMs)
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Computation of distance matrix for functional inputs
#'
#' @description Precomputes the distance matrices that will be used later inside the \emph{functional} term of the kernel function
#'              during the optimization of the hyperparmeters of the metamodel.
#'
#' @param fpIn1 a first list with as many elements as functional inputs. The i-th element must be a matrix with the projection coefficients
#' for the i-th functional input.
#' @param fpIn2 a second list with as many elements as functional inputs. The i-th element must be a matrix with the projection coefficients
#' for the i-th functional input.
#' @param J a list with as many elements as functional inputs. The i-th element must be the Gram matrix of the basis functions used for
#' the projection of the i-th functional input.
#' @return A list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#'
#' @importFrom stats cov
#' @keywords internal
#'
#' @author José Betancourt
#' @export
setFunDistance <- function(fpIn1, fpIn2, J){
  df <- length(fpIn1)
  fMs <- list()
  for (l in 1:df) {
    matrDfL1 <- as.matrix(fpIn1[[l]])
    matrDfL2 <- as.matrix(fpIn2[[l]])
    Jl <- J[[l]]
    DFbar1 <- matrDfL1 %*% Jl
    DFbar2 <- matrDfL2 %*% Jl
    Qtmp <- matrix(0L, nrow = nrow(matrDfL1), ncol = nrow(matrDfL2))
    p <- nrow(Jl)
    for (a in 1:p) {
      Qtmp <- Qtmp + (outer(matrDfL1[,a], matrDfL2[,a], '-') * outer(DFbar1[,a], DFbar2[,a], '-'))
    }
    Qtmp[Qtmp < 0] <- 0 # To prevent numerical problems like slightly negative distances
    fMs[[l]] <- sqrt(Qtmp)
  }
  return(fMs)
}
# -------------------------------------------------------------------------------------------------------------------------------------
