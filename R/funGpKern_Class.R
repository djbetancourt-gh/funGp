#' Class: data structures related to the kernel of a funGp model
#'
#' Fill this!!!!!!!!!
#'
#' @slot kerType Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#' @slot disType Object of class \code{"character"}. Distance type. To be chosen from {"scalar", "functional"}.
#' @slot varHyp Object of class \code{"numeric"}. Estimated variance parameter.
#' @slot lsHyps Object of class \code{"numeric"}. Estimated length-scale parameters.
#'
#' @rdname kernel-class
#'
#' @author José Betancourt
#' @export
setClass("funGpKern",
         representation(
                        kerType = "character",          # kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}
                        disType = "character",          # distance type. To be chosen from {"scalar", "functional"}
                        varHyp = "numeric",             # estimated variance parameter
                        lsHyps = "numeric"              # estimated length-scale parameters
                        ),
         validity = function(object) {T})
