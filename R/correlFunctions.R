#' @title Fill this!!!
#' @description Fill this!!!!!!!!!!
#'
#' @param thetas Fill this!!!!!!!!!!
#' @param Ms a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param kerType Fill this!!!!!!!!!!
#' @return Fill this!!!!!!!!!!
#'
#' @keywords internal
#'
#' @author José Betancourt
#' @export
setR <- function(thetas, Ms, kerType) {
  switch(kerType,
         "gauss" = {# 1: Gaussian
           R <- gaussian_cor(Ms, thetas)
         },

         "matern5_2" = {# 2: Matern 5/2
           R <- matern52_cor(Ms, thetas)
         },

         "matern3_2" = {# 3: Matern 3/2
           R <- matern32_cor(Ms, thetas)
         }
  )
  return(R)
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Gaussian correlation function
#' @description Applies the Gaussian correlation function to a list of distance matrices using corresponding length-scale parameters.
#'
#' @param Ms a list with as many elements as inputs. Each element of the list must be a n times n matrix of differences
#'           between the observation coordinates, with n the number of input coordinates.
#' @param thetas an array with the length-scale parameters corresponding to the inputs for which the distance matrices were provided.
#' @return A correlation matrix of dimension n times n computed from the Gaussian correlation function.
#'
#' @keywords internal
#'
#' @author José Betancourt
#' @export
gaussian_cor <- function(Ms, thetas) {
  # browser()
  sPL2 <- matrix(0L, nrow = nrow(Ms[[1]]), ncol = ncol(Ms[[1]]))
  for (l in 1:length(thetas)) {
    sPL2 <- sPL2 + (as.matrix(Ms[[l]])^2)/(thetas[l]^2)
  }
  return(exp(-0.5*sPL2))
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Matern 5/2 correlation function
#' @description Applies the Matern 5/2 correlation function to a list of distance matrices using corresponding length-scale parameters.
#'
#' @param Ms a list with as many elements as inputs. Each element of the list must be a n times n matrix of differences
#'           between the observation coordinates, with n the number of input coordinates.
#' @param thetas an array with the length-scale parameters corresponding to the inputs for which the distance matrices were provided.
#' @return A correlation matrix of dimension n times n computed from the Matern 5/2 correlation function.
#'
#' @keywords internal
#'
#' @author José Betancourt
#' @export
matern52_cor2 <- function(Ms, thetas) {
  Rm <- matrix(0L, nrow = nrow(Ms[[1]]), ncol = ncol(Ms[[1]]))
  for (l in 1:length(thetas)) {
    Rm <- Rm + (as.matrix(Ms[[l]])^2)/(thetas[l]^2)
  }
  return((1 + sqrt(5)*sqrt(Rm) + 5*Rm/3) * exp(-sqrt(5)*sqrt(Rm)))
}
# -------------------------------------------------------------------------------------------------------------------------------------

matern52_cor <- function(Ms, thetas) {
  Rm <- matrix(1L, nrow = nrow(Ms[[1]]), ncol = ncol(Ms[[1]]))
  for (l in 1:length(thetas)) {
    Rm <- Rm * (1 + sqrt(5)*abs(as.matrix(Ms[[l]]))/thetas[l] +
                  (5/3)*(as.matrix(Ms[[l]])^2)/(thetas[l]^2)) * exp(-sqrt(5)*abs(as.matrix(Ms[[l]])/thetas[l]))
  }
  return(Rm)
}

#' @title Matern 3/2 correlation function
#' @description Applies the Matern 5/2 correlation function to a list of distance matrices using corresponding length-scale parameters.
#'
#' @param Ms a list with as many elements as inputs. Each element of the list must be a n times n matrix of differences
#'           between the observation coordinates, with n the number of input coordinates.
#' @param thetas an array with the length-scale parameters corresponding to the inputs for which the distance matrices were provided.
#' @return A correlation matrix of dimension n times n computed from the Matern 3/2 correlation function.
#'
#' @keywords internal
#'
#' @author José Betancourt
#' @export
matern32_cor <- function(Ms, thetas) {
  Rm <- matrix(0L, nrow = nrow(Ms[[1]]), ncol = ncol(Ms[[1]]))
  for (l in 1:length(thetas)) {
    Rm <- Rm + (as.matrix(Ms[[l]])^2)/(thetas[l]^2)
  }
  return((1 + sqrt(3) * sqrt(Rm)) * exp(- sqrt(3) * sqrt(Rm)))
}
