# ==========================================================================================================
#' @title Printing methods for the funGp package
#' @description This set of method enables printing of the main objects defined in the funGp package. That
#'   corresponds to \linkS4class{fgpm}, \linkS4class{fgpKern}, \linkS4class{fgpProj}, and \linkS4class{Xfgpm}
#'   objects, representing funGp models, data structures related to the kernel of the model, data structures
#'   related to projection of inputs, and structures related to structural optimization of the model,
#'   respectively.
#'
#' @param object either a \linkS4class{fgpm}, \linkS4class{fgpKern}, \linkS4class{fgpProj}, or
#'   \linkS4class{Xfgpm} object.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @name show
#' @rdname show-methods
#' @importFrom methods show
if(!isGeneric("show")) {setGeneric(name = "show", def = function(object) standardGeneric("show"))}
# ==========================================================================================================
