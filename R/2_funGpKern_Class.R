# ==========================================================================================================
# Class for kernels of funGp models
# ==========================================================================================================



# ==========================================================================================================
# Developer oriented methods
# ==========================================================================================================

# Constructor of the class
# ----------------------------------------------------------------------------------------------------------
#' @title Class: data structures related to the kernel of a funGp model
#' @description Fill this!!!!!!!!!
#'
#' @slot kerType Object of class \code{"character"}. Kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}.
#' @slot disType Object of class \code{"character"}. Distance type. To be chosen from {"scalar", "functional"}.
#' @slot varHyp Object of class \code{"numeric"}. Estimated variance parameter.
#' @slot s_lsHyps Object of class \code{"numeric"}. Estimated length-scale parameters for scalar inputs.
#' @slot f_lsHyps Object of class \code{"numeric"}. Estimated length-scale parameters for functional inputs.
#'
#' @rdname kernel-class
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setClass("funGpKern",
         representation(
                        kerType = "character",          # kernel type. To be chosen from {"gauss", "matern5_2", "matern3_2"}
                        disType = "character",          # distance type. To be chosen from {"scalar", "functional"}
                        varHyp = "numeric",             # estimated variance parameter
                        s_lsHyps = "numeric",           # estimated length-scale parameters for scalar inputs
                        f_lsHyps = "numeric"            # estimated length-scale parameters for functional inputs
                        ),
         validity = function(object) {T})
# ----------------------------------------------------------------------------------------------------------


# ==========================================================================================================
# User oriented methods. For documentation of generic methods check the extraDoc.R file
# ==========================================================================================================

# Method to print the kernel of a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name show
#' @description This is my description
#' @rdname show-methods
#' @importFrom methods show
#' @param object An object to show.
#' @author José Betancourt, François Bachoc and Thierry Klein
if(!isGeneric("show")) {setGeneric(name = "show", def = function(object) standardGeneric("show"))}

#' @title Fill!!!!!!!!!!!
#' @name show
#' @rdname show-methods
#' @aliases show,funGpKern-method
# @keywords internal
setMethod("show", "funGpKern", function(object) show.funGpKern(object))

show.funGpKern <- function(object) {
  cat(paste("* Kernel type: ", object@kerType, "\n", sep = ""))
  cat(paste("* Distance type: ", object@disType, "\n\n", sep = ""))

  cat("* Hyperparameters:\n")
  cat(paste("  -> variance: ", format(object@varHyp, digits = 3, nsmall = 4), "\n", sep = ""))
  cat("  -> length-scale:\n")
  ds <- length(object@s_lsHyps)
  if (ds > 0) {
    for (i in 1:ds) {
      cat(paste("\t ls(X", i, "): ", format(object@s_lsHyps[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
  }
  df <- length(object@f_lsHyps)
  if (df > 0) {
    for (i in 1:df) {
      cat(paste("\t ls(F", i, "): ", format(object@f_lsHyps[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
  }
}
# ----------------------------------------------------------------------------------------------------------
