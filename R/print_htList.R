#' @exportS3Method print htList
print.htList <- function(x, ...) {
  cat(paste0("Heads and tails for ", ifelse(attr(x, "r"), "", "non-"), "recursive parametrization\n"))
  
  if (isTRUE(attr(x, "by_dist"))) {
    x <- purrr::transpose(x)
    x <- lapply(x, function(y) unlist(y, recursive = FALSE))
  }
  
  tmp <- rapply(x, function(x) paste(x, collapse=", "), how = "list")
  tmp <- data.frame(lapply(tmp, unlist))
  if (nrow(tmp) == 0) {
    cat("(empty graph)")
    return(invisible(x))
  }
  tmp <- data.frame(heads=tmp[,1], "."="|", tails=tmp[,2], ".."="|", intr.=tmp[,3])
  tmp$tails[nchar(tmp$tails) == 0] <- "-"
  
  cat("\n")
  print(tmp, row.names=FALSE)
  
  invisible(x)
}
