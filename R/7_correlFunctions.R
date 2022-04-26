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
         })
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
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
#' @export
gaussian_cor <- function(Ms, thetas) {
  Rm <- matrix(1L, nrow = nrow(Ms[[1]]), ncol = ncol(Ms[[1]]))
  for (l in 1:length(thetas)) {
    Rm <- Rm * exp(-0.5*(as.matrix(Ms[[l]])/thetas[l])^2)
  }
  return(Rm)
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
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
#' @export
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
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
#' @export
matern32_cor <- function(Ms, thetas) {
  Rm <- matrix(1L, nrow = nrow(Ms[[1]]), ncol = ncol(Ms[[1]]))
  for (l in 1:length(thetas)) {
    Rm <- Rm * (1 + sqrt(3)*abs(as.matrix(Ms[[l]]))/thetas[l]) * exp(-sqrt(3)*abs(as.matrix(Ms[[l]])/thetas[l]))
  }
  return(Rm)
}
