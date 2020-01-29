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
