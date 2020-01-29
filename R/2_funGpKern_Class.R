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
#' @slot f_disType Object of class \code{"character"}. Distance type. To be chosen from {"scalar", "functional"}.
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
                        f_disType = "character",        # distance type. To be chosen from {"scalar", "functional"}
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
#' @importFrom knitr kable
#' @param object An object to show.
#' @author José Betancourt, François Bachoc and Thierry Klein
if(!isGeneric("show")) {setGeneric(name = "show", def = function(object) standardGeneric("show"))}

#' @title Fill!!!!!!!!!!!
#' @name show
#' @rdname show-methods
#' @aliases show,funGpKern-method
# @keywords internal
setMethod("show", "funGpKern", function(object) show.funGpKern(kernel = object))

show.funGpKern <- function(kernel) {
  mainTxt <- "Kernel structure"
  cat(paste("\n", mainTxt, paste(rep("_", 9), collapse = ""), sep = ""))

  cat(paste("\n\n* Kernel type: ", kernel@kerType, "\n", sep = ""))
  cat("* Scalar distance: scalar\n")

  df <- length(kernel@f_lsHyps)
  if (df > 0) {
    cat("* Functional distance:")
    np <- min(df, 8)
    G <- cbind(paste("F", 1:np, sep = ""), kernel@f_disType)
    colnames(G) <- c("Input", "Distance")
    if (np < df) {
      G <- rbind(G, rep("...", 2))
    }
    print(kable(G, align = 'c', row.names = F))
  }

  cat("\n* Hyperparameters:\n")
  cat(paste("  -> variance: ", format(kernel@varHyp, digits = 3, nsmall = 4), "\n", sep = ""))
  cat("  -> length-scale:\n")
  ds <- length(kernel@s_lsHyps)
  if (ds > 0) {
    for (i in 1:ds) {
      cat(paste("\t ls(X", i, "): ", format(kernel@s_lsHyps[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
  }
  if (df > 0) {
    for (i in 1:df) {
      cat(paste("\t ls(F", i, "): ", format(kernel@f_lsHyps[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
  }
  cat(paste(rep("_", 25), collapse = ""))
}
# ----------------------------------------------------------------------------------------------------------
