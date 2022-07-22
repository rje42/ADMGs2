##' @export
intrinsicSets3 <- function (graph, r=TRUE, by_district=FALSE, sort=1, 
                            indep=FALSE, topOrd) {

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
    out <- lapply(subgrs, intrinsicSets3, r=r, by_district=FALSE, sort=sort)
    if (length(clq) > 0) out <- c(list(clq), out)
    return(out)
  }
  
  ## get a topological order
  if (missing(topOrd)) topOrd <- topologicalOrder(graph)
  
  ## get initial segment districts according to topOrd
  initDistGrs <- lapply(seq_along(subgrs), function(x) {
    subTopOrd <- intersect(topOrd, subgrs[[x]]$v)
    lapply(seq_along(subTopOrd), 
           function(i) {
             vs <- intersect(subTopOrd[seq_len(i)], dis(subgrs[[x]],subTopOrd[i]))
             tmp <- subgrs[[x]][vs]
             if (r) tmp[dis(tmp, subTopOrd[i])]
             else latentProject(tmp, latent=setdiff(graph$v, dis(tmp, subTopOrd[i])), only_directed = TRUE)
           })
  })
  initDistGrs <- unlist(initDistGrs, recursive = FALSE)

  ord <- order(sapply(initDistGrs, function(x) {
    max(match(x$v, topOrd))
  }))
  initDistGrs <- initDistGrs[ord]
  if (indep) {
    precs <- lapply(seq_along(topOrd), function(x) topOrd[seq_len(x)])
    mbs <- mapply(function(x,y) mb(graph, v=x, A=y), x=topOrd, y=precs, SIMPLIFY = FALSE)
    ind <- mapply(function(x,y) setdiff(x,y), precs, mbs, SIMPLIFY = FALSE)
    mbs <- mapply(function(x,y) setdiff(x,y), mbs, topOrd, SIMPLIFY = FALSE)
    
    inds <- mapply(as.ci, topOrd, ind, mbs, SIMPLIFY = FALSE)
  }
  
  ## now obtain intrinsic sets
  out2 <- mapply(function(x,y) intSets(x, y, topOrd, r=r, indep=indep), initDistGrs, seq_along(topOrd), SIMPLIFY = FALSE)
    # out <- lapply(subgrs, intSets, topOrd)
  # }
  
  out <- unlist(out2, recursive = FALSE)
  out <- c(clq, out)
  
  if (indep) {
    indeps <- lapply(out2, attr, "indep")
    indeps <- unlist(indeps, recursive = FALSE)
    inds <- c(inds, indeps)
    inds <- inds[sapply(inds, function(x) length(x[[2]])) > 0]
    attr(out, "indep") <- inds
  }
  
  out2 <- out

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
  if (indep) attr(out, "indep") <- attr(out2, "indep")
  
  return(out)
}

intSets <- function(graph, nt_rmv, topOrd, r=TRUE, indep=FALSE) {
  # print(graph$v)
  # cat(paste0(rep(" - ", depth)))
  # cat(graph$v)
  # cat("\n")

  ## see which vertices can be removed
  max_v <- topOrd[last(nt_rmv)]
  # if (max_v == 7) browse()
  rmv <- match(setdiff(sterile(graph[dis(graph, max_v)]), topOrd[nt_rmv]), topOrd)
  if (length(rmv) == 0) return(list(graph$v))


  # if (r) {
  dists <- lapply(seq_along(rmv), 
                  function(i) dis(graph[-topOrd[rmv[i]]], max_v, sort=2))
  if (indep) {
    pa_dists <- lapply(dists, function(x) adj(graph, x, etype="directed", dir=-1,
                                              inclusive=FALSE))
    mbs <- mapply(union, dists, pa_dists, SIMPLIFY = FALSE)
    ind <- mapply(function(x,y) setdiff(graph$v, c(x,y)), mbs, topOrd[rmv], SIMPLIFY = FALSE)
    mbs <- lapply(mbs, function(x) setdiff(x, max_v))
    inds <- mapply(as.ci, rep(max_v, length(dists)), ind, mbs, SIMPLIFY = FALSE)
  }
  dup <- duplicated(dists)
  if (any(dup)) stop()  # shouldn't be any repetition
  
  # nt_rmv_list <- rep(list(nt_rmv), length(dists))
  # dists <- dists[!dup]
  # nt_rmv_list <- mapply(function(i,x) c(topOrd[rmv[i]], nt_rmv), 
  #                                             seq_along(rmv), dists, SIMPLIFY = FALSE)
  
  
  if (r) {
    ## just take the district for recursive intrinsic sets
    nt_rmv_list <- mapply(function(i,x) match(intersect(topOrd[seq_len(length(topOrd)-rmv[i])+rmv[i]], x), topOrd),
                        seq_along(rmv), dists, SIMPLIFY = FALSE)
    out2 <- c(list(list(graph$v)), mapply(function(x,y) intSets(graph[x],y,topOrd,r=r,indep=indep), dists, nt_rmv_list, SIMPLIFY = FALSE))
  }
  else {
    ## need to include ancestors in nt_rmv for non-recursive sets
    nt_rmv_list <- mapply(function(i,x) match(intersect(topOrd[seq_len(length(topOrd)-rmv[i])+rmv[i]], anc(graph, x)), topOrd),
                        seq_along(rmv), dists, SIMPLIFY = FALSE)
    grs <- lapply(dists, function(x) latentProject(graph, v=x, only_directed = TRUE))
    out2 <- c(list(list(graph$v)), mapply(function(x,y) intSets(x,y,topOrd,r=r,indep=indep), grs, nt_rmv_list, SIMPLIFY = FALSE))
  }
  # graph2 <- graph[-graph$v[rmv]]
  # graph2 <- graph[dis(graph, max_v)]
  # out <- Recall(graph2, topOrd, nt_rmv=nt_rmv)
  
  # if (r) {
  # }
  # else {
  #   out <- c(list(list(graph$v)), mapply(function(x,y,gr) intSets(gr[x],y,topOrd,r=FALSE), x=dists, y=nt_rmv_list, gr=grs, SIMPLIFY = FALSE))
  # }
  
  out <- unlist(out2, recursive = FALSE)
  
  if (indep) {
    indeps <- lapply(out2, attr, "indep")
    indeps <- unlist(indeps, recursive = FALSE)
    attr(out, "indep") <- c(inds, indeps)
  }
  
  return(out)
}

