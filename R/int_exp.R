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
  if (r) {
    subgrs <- lapply(dists, function(x) graph[x])
  }
  else {
    subgrs <- vector(mode="list", length=length(dists))
    for (d in seq_along(dists)) {
      subgrs[[d]] <- latentProject(graph, v = dists[[d]], only_directed = TRUE)
      # for (i in dists[[d]]) {
      # 
      #   des <- MixedGraphs::pathConnected(graph[anc(graph, dists[[d]])], i, 
      #                                     setdiff(dists[[d]], i), 
      #                                     etype="directed", dir=-1)
      #   addE <- eList(lapply(des, function(x) c(i,x)))
      #   subgrs[[d]] <- addEdges(subgrs[[d]], directed=addE)
      # }
    }
  }
  
  ## if requested to give answer by district, then recurse function
  if (by_district) {
    out <- lapply(subgrs, intrinsicSets3, r=r, by_district=FALSE, sort=sort, safe=TRUE)
    if (length(clq) > 0) out <- c(list(clq), out)
    return(out)
  }
  
  ## get a topological order
  topOrd <- topologicalOrder(graph)
  
  ## get initial segment districts according to topOrd
  initDistGrs <- lapply(seq_along(subgrs), function(x) {
    subTopOrd <- intersect(topOrd, subgrs[[x]]$v)
    lapply(seq_along(subTopOrd), 
           function(i) {
             tmp <- subgrs[[x]][subTopOrd[seq_len(i)]]
             tmp[dis(tmp, subTopOrd[i])]
           })
  })
  initDistGrs <- unlist(initDistGrs, recursive = FALSE)

  ord <- order(sapply(initDistGrs, function(x) {
    max(match(x$v, topOrd))
  }))
  initDistGrs <- initDistGrs[ord]
  
  ## now obtain intrinsic sets
  out <- mapply(function(x,y) intSets(x, y, topOrd, r=r), initDistGrs, topOrd, SIMPLIFY = FALSE)
    # out <- lapply(subgrs, intSets, topOrd)
  # }
  
  out <- unlist(out, recursive = FALSE)
  out <- c(clq, out)

  ## apply sort criteria  
  if (sort > 0) {
    idx <- sapply(out, function(x) sum(2^(x-1)))
    out <- out[!duplicated(idx)]
 
    if (sort > 1) {
      out <- lapply(out, sort.int)
      
      if (sort > 2) {
        idx <- idx[!duplicated(idx)]
        out <- out[order(idx)]
      }
    }
  }
  
  return(out)
}

intSets <- function(graph, nt_rmv, topOrd, r=TRUE) {
  # print(graph$v)
  ## see which vertices can be removed
  max_v <- last(nt_rmv)
  rmv <- match(setdiff(sterile(graph[dis(graph, max_v)]), nt_rmv), topOrd)
  if (length(rmv) == 0) return(list(graph$v))


  # if (r) {
    dists <- lapply(seq_along(rmv), 
                    function(i) dis(graph[-topOrd[rmv[i]]], max_v, sort=2))
  # }
  # else {
  #   dists <- lapply(seq_along(rmv), ### USE LATENT PROJECTION
  #                   function(i) dis(graph[-topOrd[rmv[i]]], max_v, sort=2))
  # }
  dup <- duplicated(dists)
  if (any(dup)) stop()  # shouldn't be any repetition
  
  # nt_rmv_list <- rep(list(nt_rmv), length(dists))
  # dists <- dists[!dup]
  # nt_rmv_list <- mapply(function(i,x) c(topOrd[rmv[i]], nt_rmv), 
  #                                             seq_along(rmv), dists, SIMPLIFY = FALSE)
  
  
  if (r) {
    ## just take the district for recursive intrinsic sets
    nt_rmv_list <- mapply(function(i,x) intersect(topOrd[seq_len(length(topOrd)-rmv[i])+rmv[i]], x),
                        seq_along(rmv), dists, SIMPLIFY = FALSE)
    out <- c(list(list(graph$v)), mapply(function(x,y) intSets(graph[x],y,topOrd,r=TRUE), dists, nt_rmv_list, SIMPLIFY = FALSE))
  }
  else {
    ## need to include ancestors in nt_rmv for non-recursive sets
    nt_rmv_list <- mapply(function(i,x) intersect(topOrd[seq_len(length(topOrd)-rmv[i])+rmv[i]], anc(graph, x)),
                        seq_along(rmv), dists, SIMPLIFY = FALSE)
    grs <- lapply(dists, function(x) latentProject(graph, v=x, only_directed = TRUE))
    out <- c(list(list(graph$v)), mapply(function(x,y) intSets(x,y,topOrd,r=FALSE), grs, nt_rmv_list, SIMPLIFY = FALSE))
  }
  # graph2 <- graph[-graph$v[rmv]]
  # graph2 <- graph[dis(graph, max_v)]
  # out <- Recall(graph2, topOrd, nt_rmv=nt_rmv)
  
  # if (r) {
  # }
  # else {
  #   out <- c(list(list(graph$v)), mapply(function(x,y,gr) intSets(gr[x],y,topOrd,r=FALSE), x=dists, y=nt_rmv_list, gr=grs, SIMPLIFY = FALSE))
  # }
  
  return(unlist(out, recursive = FALSE))
}

