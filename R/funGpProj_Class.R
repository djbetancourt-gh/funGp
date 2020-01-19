#' Class: data structures related to projection of functional inputs
#'
#' Fill this!!!!!!!!!
#'
#' @slot doProj Object of class \code{"logical"}. Should projection of functional inputs be done?
#' @slot fpDims Object of class \code{"numeric"}. Projection dimension of each functional input.
#' @slot basis Object of class \code{"list"}. Each element (fDims_i x fpDims_i) contains the basis
#'                                            functions used for the projection of one functional input.
#' @slot coefs Object of class \code{"list"}. Each element (n x fpDims_i) contains the coefficients
#'                                            used for the projection of one functional input.
#'
#' @rdname proj-Class
#'
#' @author José Betancourt
#' @export
setClass("funGpProj",
         representation(
                        doProj = "logical",          # should projection of functional inputs be done?
                        fpDims = "numeric",          # projection dimension of each functional input
                        basis = "list",              # each element (fDims_i x fpDims_i) contains the basis
                                                     # functions used for the projection of one fun. input
                        coefs = "list"               # each element (n x fpDims_i) contains the coefficients
                                                     # used for the projection of one fun. input
                        ),
         validity = function(object) {T})
