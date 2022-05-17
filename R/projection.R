##' Project out latent variables
##' 
##' @param graph an object of class \code{mixedgraph}
##' @param latent vertices to project out
##' @param v alternatively, vertices to keep
##' @param only_directed logical: should only directed edges be added in?
##' 
##' @details \code{only_directed} should only be used in special situations,
##' such as when computing non-recursive intrisic sets.
##' 
##' @export
latentProject <- function(graph, latent, v, only_directed=FALSE, sort=1) {
  if (!is.mixedgraph(graph)) stop("'graph' must be an object of class 'mixedgraph'")
  if (!xor(missing(latent), missing(v))) stop("Must specify 'latent' or 'v' but not both")
  if (length(un(graph)) > 0) stop("Can't handle undirected edges yet")
  
  if (missing(latent)) latent <- setdiff(graph$v, v)
  else latent <- intersect(graph$v, latent)
  
  for (l in latent) {
    new_dir <- new_bi <- adjList(n=length(graph$vnames))
    
    pas <- pa(graph, l)
    chs <- ch(graph, l)
    if (only_directed) sibs <- integer(0)
    else sibs <- sib(graph, l)
    
    new_dir[chs] <- pas
    if (sort > 1) new_dir <- mapply(function(x,y) sort.int(x), new_dir)
    if (!only_directed) {
      new_bi[chs] <- list(c(sibs, chs))
      new_bi[sibs] <- list(chs)
      new_bi <- mapply(function(x,y) unique.default(setdiff(x,y)), new_bi, seq_along(new_bi))
      if (sort > 1) new_bi <- mapply(function(x,y) sort.int(x), new_bi)
    }
    class(new_dir) <- class(new_bi) <- "adjList"
    graph <- addEdges(graph, makeEdgeList(directed=new_dir, bidirected=new_bi))[-l]
  }
  
  return(graph)
}