##' Compute the (maximal) arid projection of a graph
##' 
##' @param graph an ADMG as an object of class \code{mixedgraph}
##' @param maximal should the projection computed be maximal?
##' 
##' @details Algorithm: compute the intrinsic closure of every 
##' vertex. Use this to obtain directed edges. Then add in bidirected
##' edges for vertices not already adjacent.  Then go through pairs
##' and check if intrinsic closures are joined by a bidirected edge.
##' If maximal graph asked for, then also check intrinsic closures
##' of the pairs.
##' 
aridProj <- function (graph, maximal=TRUE) {
  
  intSets <- intrinsicSets(graph)
  intSets_lens <- lengths(intSets)
  intClo <- list()
  
  ## start with existing directed edges
  dir_adj <- dir_adj_out <- adjMatrix(graph$edges$directed, directed = TRUE)
  
  ## start by getting intrinsic closure of each vertex
  for (v in graph$v) {
    if (list(v) %in% intSets) {
      ## if {v} is an intrinsic set then its own closure
      intClo[[v]] <- v
    }
    else {
      ## otherwise intrinsic closure is smallest 
      ## intrinsic set containing v
      wh <- sapply(intSets, function(x) v %in% x)
      intClo[[v]] <- intSets[wh][[which.min(intSets_lens[wh])]]
      dir_adj_out[,v] <- 1*(rowSums(dir_adj_out[,intClo[[v]],drop=FALSE]) > 0)
    }
  }
  
  ## fill in bidirected edges not replaced by directed edges
  bi_adj <- bi_adj_out <- adjMatrix(graph$edges$bidirected)
  bi_adj_out[dir_adj_out > 0 | t(dir_adj_out > 0)] = 0
 
  ## go through all possible edges to see if bidirected
  ## edge should be added...
  ## start by getting intrinsic closure of each vertex
  for (v in graph$v) for (w in graph$v[graph$v < v]) {
    if (dir_adj_out[v,w] > 0 || dir_adj_out[w,v] > 0 || bi_adj_out[v,w] > 0) {
      next
    }
    
    ## if intrinsic closures are joined by a bidirected edge,
    ## then add in
    if (any(bi_adj[intClo[[v]],intClo[[w]]] > 0)) {
      bi_adj_out[v,w] <- bi_adj_out[w,v] <- 1
      next
    }
    
    ## if we want a maximal graph, check if <v,w> is bidirected connected
    if (bi_adj_out[v,w] == 0 && maximal) {
      if (length(districts(graph[intrinsicClosure(graph, c(v,w))])) == 1) {
        bi_adj_out[v,w] <- bi_adj_out[w,v] <- 1
      }
    }
    
  }

  ## construct a graph with the new edges
  gr_out <- graph
  gr_out$edges <- list(directed=dir_adj_out, bidirected=bi_adj_out)
  return(gr_out)
}


##' Check if a graph is arid
##' 
##' @param graph ADMG of class \code{mixedgraph}
##' 
is_arid <- function(graph){
  
  w <- c(ch(graph, graph$v), sib(graph, graph$v))
  vs <- unique.default(w)
    
  ## if intrinsic closure contains more than v, then not arid
  if (length(vs) > 0) for (v in vs) {
    if (length(intrinsicClosure(graph, v)) > 1) return(FALSE)
  }
  
  return(TRUE)
}
