# ==========================================================================================================
# Extra documentation for help pages of generic user-oriented methods
# ==========================================================================================================

# Method to print a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @title Printing of a funGp model
#' @description Prints the main features of a funGp model. The list of attributes printed depends on the
#'              type of model built, which can be scalar-input, functional-input or hybrid-input.
#'
#' @param object a funGp model. Might be either a funGp_S, funGp_F or funGp_SF object depending on the type
#'               of model built.
#'
#' @name show.funGpModel
#'
#' @author José Betancourt
NULL
# ----------------------------------------------------------------------------------------------------------

# Method to get the list of projected functional inputs
# ----------------------------------------------------------------------------------------------------------
#' @title List of projected functional inputs
#' @description Provides a list with the functional inputs rebuilt from the basis functions and projection
#'              coefficients used for dimension reduction.
#'
#' @param object a funGp model. Might be either a funGp_S, funGp_F or funGp_SF object depending on the type
#'               of model built.
#'
#' @name getProjfuns
#'
#' @author José Betancourt
NULL
# ----------------------------------------------------------------------------------------------------------

# Method to get the list of Gram matrices computed from the basis functions
# ----------------------------------------------------------------------------------------------------------
#' @title List of gram matrices
#' @description Provides a list with the gram matrices built from the basis functions used to project each
#'              functional input
#'
#' @param object a funGp model. Might be either a funGp_S, funGp_F or funGp_SF object depending on the type
#'               of model built.
#'
#' @name getProjgram
#'
#' @author José Betancourt
NULL
# ----------------------------------------------------------------------------------------------------------
