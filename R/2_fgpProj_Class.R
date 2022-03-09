# ==========================================================================================================
# S4 class for structures linked to projections in a \code{fgpm} model
# ==========================================================================================================
#' @title S4 class for structures linked to projections in a \code{fgpm} model
#' @description This is the formal representation for data structures linked to projection of inputs in a
#'   Gaussian process model within the \link[=funGp-package]{funGp package}.
#'
#' @slot pdims Object of class \code{"numeric"}. Projection dimension of each input.
#' @slot basType Object of class \code{"character"}. To be chosen from {"PCA", "B-splines"}.
#' @slot basis Object of class \code{"list"}. Projection basis. For functional inputs, each element
#'   (fDims_i x fpDims_i) contains the basis functions used for the projection of one functional input.
#' @slot coefs Object of class \code{"list"}. Each element (n x fpDims_i) contains the coefficients used for
#'   the projection of one functional input.
#'
#' @rdname proj-Class
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setClass("fgpProj",
         representation(
                        pdims = "numeric",          # projection dimension of each functional input
                        basType = "character",      # family of basis functions. To be chosen from {"PCA", "B-splines"}.
                        basis = "list",             # each element (fDims_i x fpDims_i) contains the basis
                                                    # functions used for the projection of one fun. input
                        coefs = "list"              # each element (n x fpDims_i) contains the coefficients
                                                    # used for the projection of one fun. input
                        ),
         validity = function(object) {TRUE})
# ==========================================================================================================



# ==========================================================================================================
# Printer
# ==========================================================================================================
#' @rdname show-methods
#' @aliases show,fgpProj-method
setMethod("show", "fgpProj", function(object) show.fgpProj(object))

show.fgpProj <- function(object) {
  if (length(object@pdims) > 0) {
    mainTxt <- "Projection structure"
    cat(paste("\n", mainTxt, paste(rep("_", 25), collapse = ""), sep = ""))

    df <- length(object@basis)
    np <- min(df, 8)
    G <- cbind(paste("F", 1:np, sep = ""), lapply(object@basis, nrow), object@pdims, object@basType)
    colnames(G) <- c("Input", "Orig. dim", "Proj. dim", "Basis")
    if (np < df) {
      G <- rbind(G, rep("...", 4))
    }
    print(kable(G, align = 'c', row.names = FALSE))
    cat(paste(rep("_", 45), collapse = ""))

  } else {
    cat(paste("The funGp model linked to this kernel does not have functional inputs.",
           "Projection structures are not defined for it.", sep = " "))
  }
}
# ==========================================================================================================
