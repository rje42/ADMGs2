# print parameters and their values
print.mparam <-
function(x, blanks=FALSE, ...) {
  names = getMparamsNames(x, blanks=blanks)
  values = unlist(x$q)

  cat(ifelse(x$r, "Recursive", "Non-recursive"), "ADMG parametrization\n", sep=" ")

  for (i in seq_along(values)) {
    cat(names[i], " = ", values[i], "\n", sep="")
  }
}
