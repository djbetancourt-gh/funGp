# ==========================================================================================================
#' @title Summary methods for the funGp package
#' 
#' @description This set of methods produce summaries for the main
#'     \linkS4class{fgpm} and \linkS4class{Xfgpm} objects representing
#'     funGp models and structures related to structural optimization
#'     of the model, respectively.
#'
#' @param object either a \linkS4class{fgpm} or a \linkS4class{Xfgpm}
#'     object.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @name summary
#' @rdname summary-methods
#' @exportMethod summary
setGeneric(name = "summary",
           def = function(object, ...) standardGeneric("summary"))
