sample.vec <- function(x, ...){
  x[sample(length(x), ...)]
}

check.int <- function(x){
  x%%1 == 0
}

c.vecInMat_Match <- function(M, v) {
  out <- sapply(M, function(x, v) isTRUE(all.equal(x, v)), v)
  any(out)
}

r.vecInMat_Match <- function(M, v) {
  out <- apply(M, MARGIN = 1, function(x, v) isTRUE(all.equal(x, v)), v)
  which(out)
}

matInList_Match <- function(L, M) {
  out <- sapply(L, function(X, M) isTRUE(all.equal(X, M)), M)
  any(out)
}

compareNA <- function(v1,v2) {
  same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

print_2C <- function(main, cnames, gap1, gap2, c1, c2){
  h1 <- paste(main, paste(rep(" ", gap1), collapse = ""), cnames[1], sep = "")
  h2 <- paste(paste(rep(" ", gap2), collapse = ""), cnames[2], sep = "")
  cat(paste(h1, h2, "\n", sep = ""))
  for (i in 1:length(c1)) {
    s1 <- paste(paste(rep(" ", (nchar(h1)-nchar(c1[i]))), sep = "", collapse = ""), c1[i], sep = "")
    s2 <- paste(paste(rep(" ", (nchar(h2)-nchar(c2[i]))), sep = "", collapse = ""), c2[i], sep = "")
    cat(paste(s1, s2, "\n", sep = ""))
  }
}

axtags <- function(x){
  if (x > 6) {
    y <- seq_len(x)
    t <- y[(x-1) %% seq_len(x) == 0]
    l <- lapply(t, function(a) seq(1,x,a))
    s <- ((x-1)/t) + 1
    l[[which(s <= 10)[1]]]
    # l
  } else {
    seq(1, x)
  }
}


# ==========================================================================================================
# S4 class for fgpm function calls
# ==========================================================================================================
#' @title S4 class for calls to the fgpm function in funGp
#' @description User reminder of the \link[funGp]{fgpm} function call.
#'
#' @slot string Object of class \code{"character"}. User call reminder in string format.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @rdname modelCall-class
#' @export
setClass("modelCall", slots = c(string = "character"), validity = function(object) {TRUE})

if(!isGeneric("show")) {setGeneric(name = "show", def = function(object) standardGeneric("show"))}
setMethod("show", "modelCall", function(object) show.modelCall(x = object@string))
show.modelCall <- function (x, ...) {
  if (copy <- !is.null(cl <- attr(x, "class"))) {
    isNQ <- cl == "noquote"
    if (copy <- any(isNQ)) {
      ox <- x
      cl <- cl[!isNQ]
      attr(x, "class") <- if (length(cl))
        cl
    }
  }
  print(x, quote = FALSE, ...)
  invisible(if (copy) ox else x)
}
# ==========================================================================================================


# ==========================================================================================================
# S4 class for fgpm_factory function calls
# ==========================================================================================================
#' @title S4 class for fgpm_factory function calls
#' @description User reminder of the \link[funGp]{fgpm} function call.
#'
#' @slot string Object of class \code{"character"}. User call reminder in string format.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @rdname factoryCall-class
#' @export
setClass("factoryCall", slots = c(string = "character"), validity = function(object) {TRUE})

if(!isGeneric("show")) {setGeneric(name = "show", def = function(object) standardGeneric("show"))}
setMethod("show", "factoryCall", function(object) show.factoryCall(x = object@string))
show.factoryCall <- function (x, ...) {
  if (copy <- !is.null(cl <- attr(x, "class"))) {
    isNQ <- cl == "noquote"
    if (copy <- any(isNQ)) {
      ox <- x
      cl <- cl[!isNQ]
      attr(x, "class") <- if (length(cl))
        cl
    }
  }
  print(x, quote = FALSE, ...)
  invisible(if (copy) ox else x)
}

## ============================================================================
##' Abbreviate the column names and the column content of a data frame
##' describing the structure of \code{Xfgpm} object.
##'
##' @note This function is not exported. It could work with a
##'     character matrix rather than a data frame. Mind that the input
##'     names of a \code{Xfgpm} object are implicitely assumed to be
##'     "X1", "X2", ... for the scalar inputs and "F1", "F2",
##'     ... while the matrix \code{sIn} and the list \code{fIn} given
##'     at the creation time may have different names!
##' 
##' @title Put the "summary" table for a \code{Xfgpm} object in a
##'     shorter format using suitable abbreviations.
##'
##' @param structData A data frame with character columns describing
##'     the structural parameters of a \code{Xfgpm} object: state of
##'     the variables (active or not), kernel used, ...
##'
##' @return A data frame with character columns corresponding to
##'     abbreviated versions of the content of \code{data}.
##' 
##' @author Yves
##' 
formatShort <- function(structData) {
    nms <- colnames(structData)
    for (i in 1:ncol(structData)) {
        if (grepl("state", tolower(nms[i]))) {
            nms[i] <- gsub("State_", "", nms[i])
            structData[ , i] <- gsub("On", "x",  structData[ , i])
            structData[ , i] <- gsub("Off", " ",  structData[ , i])
        } else if (grepl("dist", tolower(nms[i]))) {
            nms[i] <- gsub("Distance_", "D_", nms[i])
            structData[ , i] <- gsub("L2_byindex", "idx",  structData[ , i])
            structData[ , i] <- gsub("L2_bygroup", "grp",  structData[ , i])
        } else if (grepl("basis", tolower(nms[i]))) {
            nms[i] <- gsub("Prj_basis_", "Bas_", nms[i])
            structData[ , i] <- gsub("B-splines", "Bspl",  structData[ , i])
        } else if (grepl("dim", tolower(nms[i]))) {
            nms[i] <- gsub("Prj_basis_", "bas_", nms[i])
        } else if (grepl("kernel", tolower(nms[i]))) {
            nms[i] <- "Kern"
            structData[ , i] <- gsub("matern3_2", "mat32",  structData[ , i])
            structData[ , i] <- gsub("matern5_2", "mat52",  structData[ , i])
        }
    }
    names(structData) <- nms
    structData
}
