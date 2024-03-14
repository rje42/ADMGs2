##' Project out latent variables
##' 
##' @param graph an object of class \code{mixedgraph}
##' @param latent vertices to project out
##' @param v alternatively, vertices to keep
##' @param only_directed logical: should only directed edges be added in?
##' @param sort integer:1 for unique but unsorted neighbours in the \code{adjList}, 2 for sorted before being added in.
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
  
  ## eliminate each latent variable
  for (l in latent) {
    new_dir <- new_bi <- adjList(n=length(graph$vnames))
    
    pas <- pa(graph, l)
    chs <- ch(graph, l)
    if (only_directed) sibs <- integer(0)
    else sibs <- sib(graph, l)
    
    ## orient directed edges the correct way
    if (packageVersion("MixedGraphs") < "1.0.0") new_dir[chs] <- list(pas)
    else new_dir[pas] <- list(chs)
    
    if (sort > 1) new_dir <- lapply(new_dir, sort.int)
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

##' Canonical DAG for an ADMG
##' 
##' @param graph ADMG of class \code{mixedgraph}
##' @param collect logical: should cliques of bidirected edges become a single latent?
##' @param where location to insert latent variables
##' 
##' @export
canonicalDAG <- function (graph, collect=TRUE, where=1) {
  if (where == "end") where <- nv(graph) + 1
  if (where > nv(graph) + 1) where <- nv(graph) + 1
  if (collect) {
    lts <- cliques(skeleton(graph[etype="bidirected"]))
    lts <- lts[lengths(lts) > 1]
  }
  else {
    lts <- withEdgeList(graph[etype="bidirected"])$edges
  }
  
  nl <- length(lts)
  nv <- length(graph$vnames)
  
  new_dir <- c(rep(list(NULL), nv), lts)
  class(new_dir) <- "adjList"
  if (packageVersion("MixedGraphs") < "1.0.0") new_dir <- revAdjList(new_dir)  

  ## form the graph
  out <- addNodes(graph[etype="directed"], nl)
  out <- addEdges(out, directed=new_dir)
  
  ## put latent variables in requested location
  out <- out[c(seq_len(where-1),nv + seq_len(nl), seq(from=where,to=nv, length=nv-where+1)), order=TRUE]
  # out <- out[c(nv + seq_len(nl), seq_len(nv)), order=TRUE]
  
  return(out)
}
