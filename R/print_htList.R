#' @exportS3Method print htList
print.htList <- function(x, ...) {
  cat(paste0("Heads and tails for ", ifelse(attr(x, "r"), "", "non-"), "recursive parametrization\n\n"))
  
  tmp <- rapply(x, function(x) paste(x, collapse=", "), how = "list")
  tmp <- data.frame(lapply(tmp, unlist))
  tmp <- data.frame(heads=tmp[,1], V="|", tails=tmp[,2], T="|", intr.=tmp[,3])
  
  print(tmp, row.names=FALSE)
  
  invisible(x)
}
