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
  # This function returns TRUE wherever elements are the same, including NA's, and false everywhere else.
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

# ==========================================================================================================
setClass("modelCall", slots = c(string = "character"), validity = function(object) {T})

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
setClass("factoryCall", slots = c(string = "character"), validity = function(object) {T})

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
# ==========================================================================================================
