# ==========================================================================================================
# Class for projs of funGp models
# ==========================================================================================================



# ==========================================================================================================
# Developer oriented methods
# ==========================================================================================================

# Constructor of the class
# ----------------------------------------------------------------------------------------------------------
#' @title Class: data structures related to projection of functional inputs
#' @description Fill this!!!!!!!!!
#'
#' @slot pdims Object of class \code{"numeric"}. Projection dimension of each input.
#' @slot basType Object of class \code{"character"}. To be chosen from {"PCA", "B-splines"}.
#' @slot basis Object of class \code{"list"}. Projection basis. For functioanl inputs, eEach element (fDims_i x fpDims_i)
#'                                            contains the basis functions used for the projection of one functional input.
#' @slot coefs Object of class \code{"list"}. Each element (n x fpDims_i) contains the coefficients
#'                                            used for the projection of one functional input.
#'
#' @rdname proj-Class
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setClass("funGpProj",
         representation(
                        pdims = "numeric",          # projection dimension of each functional input
                        basType = "character",       # family of basis functions. To be chosen from {"PCA", "B-splines"}.
                        basis = "list",             # each element (fDims_i x fpDims_i) contains the basis
                                                    # functions used for the projection of one fun. input
                        coefs = "list"              # each element (n x fpDims_i) contains the coefficients
                                                    # used for the projection of one fun. input
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
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @param object An object to show.
if(!isGeneric("show")) {setGeneric(name = "show", def = function(object) standardGeneric("show"))}

#' @title Fill!!!!!!!!!!!
#' @name show
#' @rdname show-methods
#' @aliases show,funGpProj-method
#' @keywords internal
setMethod("show", "funGpProj", function(object) show.funGpProj(object))

show.funGpProj <- function(object) {
  if (length(object@doProj) == 1) {
    if (object@doProj) {
      cat(paste("* Do projection: ", object@doProj, "\n", sep = ""))
      if (object@doProj) {
        cat("  -> Proj. dimension:\n")
        for (i in 1:length(object@fpDims)) {
          if (object@fpDims[i] > 0) {
            cat(paste("\t F", i, ": ", object@fpDims[i], "\n", sep = ""))
          } else {
            cat(paste("\t F", i, ": not required\n", sep = ""))
          }
        }
      }
    } else {
      cat(paste("* Do projection: ", object@doProj, ". Projection was not required for this model.", sep = ""))
    }
  } else {
    cat(paste("The funGp model linked to this kernel does not have functional inputs.",
           "Projection structures are not defined for it.", sep = " "))
  }
}
# ----------------------------------------------------------------------------------------------------------
