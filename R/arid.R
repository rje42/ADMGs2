##' Compute the (maximal) arid or ancestral projection of a graph
##' 
##' @param graph an ADMG as an object of class \code{mixedgraph}
##' @param maximal should the projection computed be maximal?
##' @param verbose logical: should additional output be given?
##' 
##' @details The input graph should be a summary graph for \code{ancProj},
##' and an ADMG for \code{aridProj}. 
##' 
##' Algorithm for \code{aridProj()}: compute the intrinsic closure of every 
##' vertex. Use this to obtain directed edges. Then add in bidirected
##' edges for vertices not already adjacent.  Then go through pairs
##' and check if intrinsic closures are joined by a bidirected edge.
##' 
##' Algorithm for \code{ancProj()}: go through in a topological order and check
##' no siblings are among ancestors.  If there are at say \code{v}, then look for parents and 
##' siblings, recursively if the latter are also ancestors.  Remove any bidirected
##' edges from discovered vertices to \code{v} and add directed edges from them
##' to \code{v}
##' 
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
  else dir_adj_out <- adjMatrix(n=n)
  
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

##' @describeIn aridProj obtain ancestral projection
##' @param directed logical: is graph an ADMG?
##' @export
ancProj <- function (graph, maximal=TRUE, directed=FALSE) {
  if (nv(graph) <= 1) return(graph)
  if (!directed && !is_SG(graph)) stop("'graph' must be a summary graph of class 'mixedgraph'")
  else if (directed && !is_ADMG(graph)) stop("'graph' must be an ADMG graph of class 'mixedgraph'")
  is_anc <- is_ancestral(graph, .top_ord = TRUE, .anc = TRUE)

  ## extract side information from is_ancestral()
  top_ord <- attr(is_anc, "top_ord")
  ancs <- attr(is_anc, "ancs")
  
  while (!is_anc) {
    ## find vertex that is not ancestral
    v <- attr(is_anc, "which")
    # wh <- which.min(c(lengths(ancs[top_ord]))
    # v <- top_ord[wh]
    gr <- graph[top_ord[seq_len(match(v, top_ord))]]
    ## get its sibling-ancestors
    viol <- intersect(sib(gr, v), ancs[[v]])
    pas <- pa(gr, viol)
    sbs <- setdiff(sib(gr, viol), c(viol,v))

    ## set up loop to obtain all vertices with edges into v
    sbs_a <- new <- intersect(sbs, ancs[[v]])  ## get cases where siblings are also ancestors
    sbs_n <- setdiff(sbs, sbs_a)
    all_a <- c(viol, sbs_a, pas)
    all_n <- sbs_n
    
    while (length(new) > 0) {
      # new <- setdiff(sib(gr, sbs_a), all)  ## get cases where siblings are also ancestors
      new_sbs <- setdiff(sib(gr, sbs_a), c(all_a, all_n))
      new <- intersect(new_sbs, ancs[[v]])
      all_a <- c(all_a, new)
      all_n <- c(all_n, setdiff(new_sbs, new))
      new_pa <- setdiff(pa(gr, sbs_a), c(all_a, all_n))
      all_a <- c(all_a, new_pa)
    }
    
    ## now add in new edges and remove old ones
    graph <- addEdges(graph, makeEdgeList(directed=eList(lapply(all_a, function(x) c(x,v)))))
    graph <- removeEdges(graph, makeEdgeList(bidirected=eList(lapply(all_a, function(x) c(x,v)))), force = TRUE)
    graph <- addEdges(graph, makeEdgeList(bidirected=eList(lapply(all_n, function(x) c(x,v)))))
    
    ## check if graph is now ancestral
    is_anc <- is_ancestral(graph, top_ord=top_ord, .anc=TRUE)
    ancs <- attr(is_anc, "ancs")
    
    ## add in directed edges viol and pas -> v
    ## add in directed edges for sbs -> v if an ancestor (does this necessitate recursive approach?)
    ## add in bidirected edges for sbs <-> v if not an ancestor
    # addEdges(graph, )
  }

  if (maximal && nedge(graph, "bidirected") > 0) {
    dists <- districts(graph)
    len_d <- lengths(dists)
    to_add <- list()
    for (d in which(len_d > 2)) {
      combn(dists[[d]], 2, simplify = FALSE)
      
      for (i in seq_along(dists[[d]][-1])) for (j in seq_len(i-1)) { 
        i1 <- dists[[d]][i]; j1 <- dists[[d]][j]
        if (!(i1 %in% adj(graph, j1))) {
          if (j1 %in% dis(graph[union(ancs[[i1]], ancs[[j1]])], i1)) {
            to_add <- c(to_add, list(c(i1,j1)))
          }
        }
      }
    }
    graph <- addEdges(graph, bidirected=eList(to_add))
  }
  
  return(graph)
}

##' @describeIn is_MArG check if a graph is arid
##' @export
is_arid <- function(graph) {
  
  if (!is_SG(graph)) return(FALSE)
  vs <- intersect(ch(graph, graph$v), sib(graph, graph$v))

  ## if intrinsic closure contains more than v, then not arid
  if (length(vs) > 0) for (v in vs) {
    if (length(intrinsicClosure(graph, v)) > 1) return(FALSE)
  }
  
  return(TRUE)
}

##' @describeIn is_MArG check if graph is maximal
##' @export
is_maximal <- function(graph, check=TRUE, ancestral, failure=FALSE) {
  
  if (check) {
    if (!is_SG(graph)) {
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
        oth_ver <- setdiff(v[seq_len(i-1)], adj(graph, v[i]))
        for (j in seq_along(oth_ver)) {
          ## look at pairs of non-adjacent vertices
          con_comp <- dis(graph[anc(graph, c(v[i],oth_ver[j]))], v[i])
          if (oth_ver[j] %in% con_comp) return(FALSE)
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
          if (j %in% dis(graph[intClo], i)) {
            out <- FALSE
            if (failure) attr(out, "failure") <- c(i,j)
            return(out)
          }
        }
      }
    }
  }
  
  return(TRUE)
}


##' Check if a graph is maximal and arid
##' 
##' @param graph summary graph or ADMG of class \code{mixedgraph}
##' @param directed logical: if \code{TRUE}, undirected edges are not allowed
##' @param failure logical: if graph is not maximal, should missing edge be returned?
##' 
##' @details Checks if the graph is both maximal and arid
##' using the functions \code{is_arid} and \code{is_maximal}.
##' 
##' @export
is_MArG <- function(graph, directed=FALSE, failure=FALSE) {
  if (directed && nedge(graph, "undirected") > 0) return(FALSE)
  else if (nedge(graph, setdiff(names(graph$edges), c("undirected", "directed", "bidirected")) > 0)) return(FALSE)
  
  if(!is_arid(graph)) return(FALSE)
  out <- is_maximal(graph, check=FALSE, failure=failure)
  
  return(out)
}

##' Check if a directed mixed graph is ancestral
##' 
##' @param graph object of class \code{mixedgraph}
##' @param top_ord optionally, a topological order
##' @param .top_ord logical: should topological order be returned if computed?
##' 
##' @details \code{graph} should only contain directed and
##' bidirected edges.
##' 
##' @export
is_ancestral <- function(graph, top_ord, .top_ord=FALSE, directed=FALSE, .anc=FALSE) {
  
  ## check no arrows point to an undirected edge
  if (!directed) {
    un_g <- un(graph)
    if (length(un_g) > 0) {
      check_un <- any(ch(graph, graph$v) %in% un_g) || any(sib(graph, graph$v) %in% un_g)
      if (check_un) {
        out <- FALSE
        if (.top_ord) {
          attr(out, "top_ord") <- NA
          attr(out, "top_ord_code") <- 2
        }
        return(out)
      }
    }
  }
  ## if graph is cyclic, return FALSE, otherwise get topological order
  if (missing(top_ord)) vs <- topologicalOrder(graph, warn=FALSE)
  else vs <- top_ord
  
  if (is.na(vs[1])) {
    out <- FALSE
    if (.top_ord) {
      attr(out, "top_ord") <- NA
      attr(out, "top_ord_code") <- 1
    }
    if (.anc) {
      attr(out, "ancs") = NA
    }
    return(out)
  }

  ## get output
  out <- TRUE
  if (.top_ord) {
    attr(out, "top_ord") <- vs
    attr(out, "top_ord_code") <- 0
  }
  
  ## now check for siblings amongst ancestors
  ancs <- vector(mode="list", length=length(vnames(graph)))
  
  for (v in vs) {
    ancs[[v]] <- c(v, unlist(ancs[pa(graph, v)]))
    if (.anc) attr(out, "ancs") = ancs
    
    if (any(sib(graph, v) %in% ancs[[v]])) {
      attr(out, "which") = v
      out[] <- FALSE
      return(out)
    }
  }
  
  return(out)
}


##' Check if a graph is maximal and ancestral
##' 
##' @param graph an object of class \code{mixedgraph}
##' @param directed logical: should this be a directed MAG?
##' @param .top_ord logical: if computed by ancestral should topological order be returned?
##' 
##' @details Checks edge types, and then runs \code{is_ancestral} and 
##' \code{is_maximal}, reports the success of both or the failure of either.
##' 
##' @export
is_MAG <- function(graph, directed=FALSE, .top_ord=FALSE) {
  if (directed && nedge(graph, "undirected") > 0) return(FALSE)
  else if (nedge(graph, setdiff(names(graph$edges), c("undirected", "directed", "bidirected")) > 0)) return(FALSE)
  
  ancestral <- is_ancestral(graph, .top_ord=.top_ord)
  if (!ancestral) return(ancestral)
  
  out <- is_maximal(graph, check=FALSE, ancestral=ancestral)
  if (.top_ord) attributes(out) <- attributes(ancestral)
  
  return(out)
}
