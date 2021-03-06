##' Compute the (maximal) arid projection of a graph
##' 
##' @param graph an ADMG as an object of class \code{mixedgraph}
##' @param maximal should the projection computed be maximal?
##' @param verbose logical: should additional output be given?
##' 
##' @details Algorithm: compute the intrinsic closure of every 
##' vertex. Use this to obtain directed edges. Then add in bidirected
##' edges for vertices not already adjacent.  Then go through pairs
##' and check if intrinsic closures are joined by a bidirected edge.
##' If maximal graph asked for, then also check intrinsic closures
##' of the pairs.
##' 
##' @export
aridProj <- function (graph, maximal=TRUE, verbose=FALSE) {
  
  n <- length(graph$vnames)
  # if (verbose) cat("Getting intrinsicSets()...")
  # # intSets <- intrinsicSets(graph)
  # if (verbose) cat("done\n")
  # intSets_lens <- lengths(intSets)
  intClo <- vector(mode="list", length=n)
  intClo[graph$v] <- lapply(graph$v, function(x) intrinsicClosure(graph, x))
  
  gr_out <- graph[etype="directed"]
  
  if (nedge(graph, "directed") > 0) {
    ## start with existing directed edges
    dir_adj <- dir_adj_out <- adjMatrix(graph$edges$directed, n=n, directed = TRUE)
    
    ## start by getting intrinsic closure of each vertex
    for (v in graph$v) {
      # if (length(intClo[[v]]) == 1) next
      
      Pa_v <- pa(graph, intClo[[v]])
      if (length(Pa_v) == 0) next
      
      dir_v <- lapply(Pa_v, function(x) c(x,v))
      class(dir_v) <- "eList"
      gr_out <- addEdges(gr_out, makeEdgeList(directed=dir_v))
    }
    # if (length(intClo[v]) > 1) {
    #   if (verbose) cat(v)
    #   
    #   dir_adj_out[intClo[[v]],v] <- 1L
    #   dir_adj_out[v,v] <- 0L
    #   
    #   # if (list(v) %in% intSets) {
    #   #   ## if {v} is an intrinsic set then its own closure
    #   #   intClo[[v]] <- v
    #   # }
    #   # else {
    #   #   ## otherwise intrinsic closure is smallest 
    #   #   ## intrinsic set containing v
    #   #   wh <- sapply(intSets, function(x) v %in% x)
    #   #   intClo[[v]] <- intSets[wh][[which.min(intSets_lens[wh])]]
    #   #   dir_adj_out[,v] <- 1*(rowSums(dir_adj_out[,intClo[[v]],drop=FALSE]) > 0)
    #   # }
    #   if (verbose) cat(" ")
    # }
    # if (verbose) cat("\n")
    
    dir_adj_out <- withAdjMatrix(gr_out[etype="directed"])$edges$directed
  }
  else dir_adj_out <- adjMatrix(n)
  
  if (nedge(graph, "bidirected") > 0) {
    ## fill in bidirected edges not replaced by directed edges
    bi_adj <- bi_adj_out <- adjMatrix(graph$edges$bidirected, n=n)
    bi_adj_out[dir_adj_out > 0 | t(dir_adj_out > 0)] = 0
    
    if (maximal) dists <- districts(graph)
    
    ## go through all possible edges to see if bidirected
    ## edge should be added...
    ## start by getting intrinsic closure of each vertex
    if (verbose) cat("Looking at: ")
    for (v in graph$v) for (w in graph$v[graph$v < v]) {
      if (verbose) cat(paste0("(", v, ",", w, ")\n"))
      
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
      if (maximal) {
        if (subsetmatch(list(c(v,w)), dists, nomatch = 0) > 0) {
          if (bi_adj_out[v,w] == 0) {
            if (length(districts(graph[intrinsicClosure(graph, c(v,w))])) == 1) {
              bi_adj_out[v,w] <- bi_adj_out[w,v] <- 1
            }
          }
        }
      }
    }
  }
  else bi_adj_out <- adjMatrix(n=n)
    
  ## construct a graph with the new edges
  gr_out <- graph
  gr_out$edges <- makeEdgeList(directed=dir_adj_out, bidirected=bi_adj_out)
  return(gr_out)
}


##' @describeIn is_MArG check if a graph is arid
##' @export
is_arid <- function(graph){
  
  if (!is.SG(graph)) return(FALSE)
  vs <- intersect(ch(graph, graph$v), sib(graph, graph$v))

  ## if intrinsic closure contains more than v, then not arid
  if (length(vs) > 0) for (v in vs) {
    if (length(intrinsicClosure(graph, v)) > 1) return(FALSE)
  }
  
  return(TRUE)
}

##' @describeIn is_MArG check if graph is maximal
##' @export
is_maximal <- function(graph, check=TRUE, ancestral) {
  
  if (check) {
    if (!is.SG(graph)) {
      return(FALSE)
    }
    if (!is_arid(graph)) return(FALSE)
  }
  
  if (missing(ancestral)) ancestral <- is_ancestral(graph)
    
  ## go through districts looking for missing bidirected edges
  dis <- districts(graph)
  dis <- dis[lengths(dis) > 2]
  
  if (ancestral) {
    for (d in seq_along(dis)) {
      v <- dis[[d]]
      for (i in seq_along(v)[-1]) {
        for (j in setdiff(v[seq_len(i-1)], adj(graph, v[i]))) {
          ## look at pairs of non-adjacent vertices
          con_comp <- dis(graph[anc(graph, v[c(i,j)])], v[i])
          if (v[j] %in% con_comp) return(FALSE)
        }
      }
    }
  }
  else {
    for (d in seq_along(dis)) {
      v <- dis[[d]]
      k <- 1
      for (i in v[-1]) {
        k <- k+1
        for (j in setdiff(v[seq_len(k-1)], adj(graph, i))) {
          ## look at pairs of non-adjacent vertices
          intClo <- intrinsicClosure(graph, c(i, j))
          if (j %in% dis(graph[intClo], i)) return(FALSE)
        }
      }
    }
  }
  
  return(TRUE)
}


##' Check if a graph is maximal and arid
##' 
##' @param graph summary graph or ADMG of class \code{mixedgraph}
##' 
##' @details Checks if the graph is both maximal and arid
##' using the functions \code{is_arid} and \code{is_maximal}.
##' 
##' @export
is_MArG <- function(graph) {
  
  if(!is_arid(graph)) return(FALSE)
  if(!is_maximal(graph, check=FALSE)) return(FALSE)

  return(TRUE)
}

##' Check if a directed mixed graph is ancestral
##' 
##' @param graph object of class \code{mixedgraph}
##' 
##' @details \code{graph} should only contain directed and
##' bidirected edges.
##' 
##' @export
is_ancestral <- function(graph) {
  
  ## check no arrows point to an undirected edge
  un_g <- un(graph)
  if (length(un_g) > 0) {
    check_un <- any(ch(graph, graph$v) %in% un_g) || any(sib(graph, graph$v) %in% un_g)
    if (check_un) return(FALSE)
  }
  
  ## if graph is cyclic, return FALSE, otherwise get topological order
  vs <- topologicalOrder(graph, warn=FALSE)
  if (is.na(vs[1])) return(FALSE)

  
  ## now check for siblings amongst ancestors
  ancs <- vector(mode="list", length=length(vnames(graph)))
  
  for (v in vs) {
    ancs[[v]] <- c(v, unlist(ancs[pa(graph, v)]))
    if (any(sib(graph, v) %in% ancs[[v]])) return(FALSE)
  }
  
  return(TRUE)
}


##' Check if a graph is maximal and ancestral
##' 
##' @param graph an object of class \code{mixedgraph}
##' 
##' @details Runs \code{is_ancestral} and then 
##' \code{is_maximal}, reports the success of both 
##' or the failure of either.
##' 
##' @export
is_MAG <- function(graph) {
  ancestral <- is_ancestral(graph)
  if (!ancestral) return(FALSE)
  else return(is_maximal(graph, check=FALSE, ancestral=ancestral))
}
