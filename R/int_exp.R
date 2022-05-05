##' @export
intrinsicSets3 <- function (graph, r=TRUE, by_district=FALSE, sort=1, safe=FALSE) {

  ## check graph is a summary graph and remove undirected part
  if(!is_SG(graph)) stop("Graph appears not to be a summary graph or ADMG")
    
  un_g <- un(graph)
  if (length(un_g) > 0) {
    clq <- cliques(graph[un_g])
    graph <- graph[-un_g]
  }
  else clq <- vector(mode="list", length=0L)

  ## now look at districts of what remains
  dists <- districts(graph)
  subgrs <- lapply(dists, function(x) graph[x])
  if (!r) {
    for (d in seq_along(dists)) {
      for (i in dists[[d]]) {
        des <- MixedGraphs::pathConnected(graph[anc(graph, dists[[d]])], i, 
                                          setdiff(dists[[d]], i), 
                                          etype="directed", dir=-1)
        addE <- eList(lapply(des, function(x) c(i,x)))
        subgrs[[d]] <- addEdges(subgrs[[d]], directed=addE)
      }
    }
  }
  
  ## if requested to give answer by district, then recurse function
  if (by_district) {
    out <- lapply(subgrs, intrinsicSets3, r=r, by_district=FALSE, sort=sort, safe=TRUE)
    if (length(clq) > 0) out <- c(list(clq), out)
    return(out)
  }
  
  topOrd <- topologicalOrder(graph)
  ## get initial segment districts
  
  initDistGrs <- lapply(seq_along(subgrs), function(x) {
    subTopOrd <- intersect(topOrd, subgrs[[x]]$v)
    lapply(seq_along(subTopOrd), 
           function(i) subgrs[[x]][subTopOrd[seq_len(i)]])
  })
  initDistGrs <- unlist(initDistGrs, recursive = FALSE)

  ord <- order(sapply(initDistGrs, function(x) {
    max(match(x$v, topOrd))
  }))
  initDistGrs <- initDistGrs[ord]
  
  # dists <- districts(graph)
  # int <- dists
  
  # if (r) {
  #   # out <- lapply(dists, function(x) intSets(graph[x], topOrd))
  #   out <- mapply(function(x,y) intSets(graph[x], y, topOrd), dists, topOrd, SIMPLIFY = FALSE)
  # }
  # else {
  #   subgrs <- lapply(dists, function(x) graph[x])
  #   for (d in seq_along(dists)) {
  #     for (i in dists[[d]]) {
  #       des <- MixedGraphs:::find_and_stop(graph[anc(graph, dists[[d]])], i, 
  #                                          setdiff(dists[[d]], i), 
  #                                          etype="directed", dir=-1)
  #       addE <- eList(lapply(des, function(x) c(i,x)))
  #       subgrs[[d]] <- addEdges(subgrs[[d]], directed=addE)
  #     }
  #   }
    
  ## now obtain intrinsic sets
  out <- mapply(function(x,y) intSets(x, y, topOrd), initDistGrs, topOrd, SIMPLIFY = FALSE)
    # out <- lapply(subgrs, intSets, topOrd)
  # }
  
  out <- unlist(out, recursive = FALSE)
  out <- c(clq, out)
  
  if (sort > 1) {
    out <- lapply(out, sort.int)
    if (sort > 2) {
      ord <- order(sapply(out, function(x) sum(2^(x-1))))
      out <- out[ord]
    }
  }
  
  return(out)
}

intSets <- function(graph, nt_rmv, topOrd) {
  ## see which vertices can be removed
  rmv <- match(setdiff(sterile(graph), nt_rmv), topOrd)
  if (length(rmv) == 0) return(list(graph$v))
  max_v <- last(nt_rmv)
  
  dists <- lapply(seq_along(rmv), 
                  function(i) dis(graph[-topOrd[rmv[seq_len(i)]]], max_v))
  nt_rmv_list <- mapply(function(i,x) intersect(topOrd[seq_len(length(topOrd)-rmv[i])+rmv[i]], x), 
                        seq_along(rmv), dists)
  # graph2 <- graph[-graph$v[rmv]]
  # graph2 <- graph[dis(graph, max_v)]
  # out <- Recall(graph2, topOrd, nt_rmv=nt_rmv)
  
  out <- c(list(list(graph$v)), mapply(function(x,y) intSets(graph[x],y,topOrd), dists, nt_rmv_list, SIMPLIFY = FALSE))
  
  return(unlist(out, recursive = FALSE))
}
