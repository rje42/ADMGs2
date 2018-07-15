#' Find best ADMG model to data.
#' 
#' Fit discrete data to the Markov structure implied by an acyclic directed
#' mixed graph.
#' 
#' Given the current graph, the algorithm fits all graphs obtained by removing,
#' adding or (in the case of the full search) replacing a single edge, and
#' picks the graph which is optimal with respect to the chosen criterion (AIC
#' or BIC).  The algorithm terminates when no move results in an improvement.
#' 
#' The default starting graph has no edges.
#' 
#' @param dat Either an object of class \code{freq.table}, or a vector of
#' counts in reverse lexicographical order.
#' @param dims The number of categories for each variable.
#' @param start Graph to start search procedure at; defaults to graph with no
#' edges.
#' @param criterion Criterion to minimize: either \code{"AIC"} or \code{"BIC"}.
#' @param quietly Logical indicating whether to suppress output.
#' @param r Logical indicating whether recursive factorizations should be used.
#' @param search \code{"full"} or \code{"basic"} search.
#' @param tol tolerance used to assess convergence in \code{\link{fitADMG}}.
#' @return An object of class \code{mixed_fit}.
#' @note Note that this procedure is not guaranteed to find the best model fit
#' for all ADMGs, but will find a 'local maximum' with respect to edge addition
#' and removal.
#' @author Ilya Shpitser, Robin Evans
#' @seealso \code{\link{fitADMG}}.
#' @references Evans and Richardson (2009) - Maximum likelihood fitting of
#' acyclic directed mixed graphs to binary data.
#' @keywords graphs
#' @examples
#' 
#' data(gss.small)
#' \dontrun{out = autoFit(gss.small)}
#' 
#' \dontrun{summary(out$graph)}
#' 
autoFit <-
function(dat, start=NULL, criterion = "AIC", quietly = FALSE, r = TRUE, search = "basic", tol = 1e-2) {
  dims <- dim(dat)
  if (length(dims) > 2) {
    dat = cbind(combinations(dims), c(dat))
    n <- length(dims)
    if (!is.null(dimnames(dat))) colnames(dat) <- c(dimnames(dat), "freq")
    else colnames(dat) <- c(paste("x", seq_len(n), sep=""), "freq")
    dat <- as.data.frame(dat)
  }
  else if (length(dims) == 2) {
  # if (length(dim(dat)) > 0) n = dim(dat)[2]-1
  # else if (!missing(dims)) n = length(dims)
  # else stop("Must provde dimensions of model")
    n <- ncol(dat)-1
    dims = apply(dat, 2, function(x) max(x)+1)[seq_len(n)]
  }
  # else if (is.vector(dat)) {
  #   dat = cbind(combinations(dims), dat)
  #   vnames = character(n)
  #   for (i in seq_len(n)) vnames[i] = paste("x", i, sep="")
  #   colnames(dat)[seq_len(n)] = vnames
  # }

  # vnames = colnames(dat)[-(n+1)]

  if (is.null(start)) {
    graph = mixedgraph(n)
  }
  else graph = start

  best_fit = fitADMG(dat, graph, tol=tol, quietly=TRUE, sparse=TRUE, r = r)
  stat = switch(criterion, AIC=summary(best_fit, fisher=FALSE)$AIC, BIC=summary(best_fit, fisher=FALSE)$BIC)

#  tmp = searchADMG(emptygraph, gss, summary(f)$AIC)

  ok = TRUE

  while(ok) {

  tmp = switch(search,  
                basic = .searchADMG(best_fit$graph, dat, stat, criterion = criterion, quietly=quietly, tol=tol),
                full = .searchADMG2(best_fit$graph, dat, stat, criterion = criterion, quietly=quietly, tol=tol),
                exhaustive = stop("not implemented"))

    ok = tmp$moved
    if (ok) {
      best_fit = tmp
      stat = tmp$stat
      if (!quietly) print(best_fit)
    }
  }

  return(best_fit)
}
