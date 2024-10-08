# getSl3 <- function(graph, top=FALSE) {
#   A <- graph
#   
#   top <- topologicalOrder(graph)
#   
#   # if (!isTopological(A, A$v)) stop("algorithm won't work")
#   
#   # S <- set()
#   S <- Sh1 <- Sh2 <- list()
#   Lv <- length(A$v)
#   
#   ## deal with smaller graphs
#   if (Lv == 0) return(S)
#   if (Lv == 1) return(list(A$v))
#   
#   if (length(A$edges$directed) == 0){
#     A$edges$directed <- matrix(0,Lv,Lv)
#   }
#   if (length(A$edges$bidirected) == 0){
#     A$edges$bidirected <- matrix(0,Lv,Lv)
#   }
#   # if graphs are DAGs or contains only bidirected edges, graphCr 
#   # only gives one edge matrix so above two 'if' are to aviod such problems.
#   
#   ## time saving: get all parents before we start
#   pa <- lapply(A$v, function(i) pa(A, i))
#   anc <- list()
#   
#   ## technically the vertex is the number in A$v (and these might not be consecutive integers)
#   for (i in top) {
#     S <- c(S, list(A$v[i]))
#     #  S <- S|set(set(A$vnames[i]))
#     
#     if (length(pa[[i]]) > 0){
#       anc[[i]] <- unique.default(c(unlist(anc[pa[[i]]]), i))
#       
#       for (j in seq_along(pa[[i]])){
#         # S <- S|set(set(A$vnames[i],A$vnames[pa_i[j]]))
#         Sh1 <- c(Sh1, list(A$v[c(pa[[i]][j], i)])) 
#       }
#     }
#     else anc[[i]] <- c(i)
#     
#     if (length(pa[[i]]) > 1) {
#       for (j in seq_along(pa[[i]])[-length(pa[[i]])]){
#         for (k in seq_len(length(pa[[i]])-j) + j){
#           # S <- S|set(set(A$vnames[i],A$vnames[pa_i[j]],A$vnames[pa_i[k]]))
#           S <- c(S, list(A$v[c(pa[[i]][j], pa[[i]][k], i)])) # 
#         }
#       }
#     }
#     #Above "for" loop is to identify heads of single vertex and corresponding tail(parents). Second if is for multiple parents
#   }
#   
#   ## here we look for heads of size two and their tails
#   ## this seems to be O(|V||E|+|E|^2)
#   for (i in A$v) {
#     for (j in which(A$edges$bidirected[seq(from=i+1,to=Lv,length.out = Lv-i),i] == 1)+i) {
#       Sh2 <- c(Sh2, list(A$v[c(i,j)]))
#       
#       B <- A[c(anc[[i]], anc[[j]])]
#       dis_i <- dis(B, i)  # |V|+|E|
#       tail_ij <- grp(B, dis_i, etype = "directed", dir=-1, inclusive = TRUE) # |V|
#       tail_ij <- setdiff(tail_ij, c(i,j))
#       
#       ## add in variables from tail
#       for (k in tail_ij) S <- c(S, list(A$v[c(k,i,j)]))
#     }
#   }
#   
#   ## now we look for heads of size 3
#   for (i in A$v) {
#     wh <- which(map_lgl(Sh2, function(x) (i %in% x)))
#     # wh <- which(map_lgl(Sh2, function(x) x[2]==i))   # assuming i is last is a bad idea!
#     
#     # tmp1 <- sapply(Sh1[wh], function(x) x[1])
#     tmp <- sapply(Sh2[wh], function(x) setdiff(x, i))
#     
#     if (length(tmp) >= 2) { 
#       ## if 2 out of 3 are heads, then can add in c(i,j,k)
#       tmp2 <- matrix(tmp[combn(length(tmp), 2)], nrow=2)
#       for (j in seq_len(ncol(tmp2))) S <- c(S, list(c(i,tmp2[,j])))
#     }
#     
#     ## if only 1 out of 3 is a head, then need to check 
#     ## if c(i,j,d) is a head for each d in dis(i)
#     for (w in seq_along(wh)) {
#       j <- Sh2[[wh[w]]][1]
#       wh2 <- which(map_lgl(Sh2, function(x) i %in% x || j %in% x)) 
#       wh2 <- setdiff(dis(A,i), c(unlist(Sh2[wh2]), i, j))
#       
#       # j <- setdiff(unlist(Sh2[wh]), i)
#       # B = subGraph(A, c(anc[[i]], anc[[j]]))
#       # disij <- dis(B, i)
#       # tl <- setdiff(grp(B, v = disij, etype="directed", dir=-1, inclusive = TRUE), c(i,j))
#       # S <- c(S, lapply(tl, function(x) c(x,i,j)))
#       
#       #      dis_i <- dis(A, i)
#       for (d in wh2) {
#         B <- A[c(anc[[i]], anc[[j]], anc[[d]])]  # |V|+|E|
#         if (d %in% dis(B, i)) S <- c(S, list(c(i,j,d)))  
#       }
#     }
#   }
#   
#   S <- c(S, Sh1, Sh2)
#   
#   #   if (i < Lv) {
#   #     for (j in seq_len(Lv-i)+i) {
#   #       if (A$edges$bidirected[i,j] == 1) {
#   #         S <- c(S, list(A$v[c(i,j)]))
#   #       }
#   #     }
#   #   }
#   # }
#   # 
#   #          
#   #         B = subGraph(A,c(anc[[i]], anc[[j]]))
#   #         disij <- dis(B, i)
#   #         distail <- disij[(disij != i) & (disij != j)]
#   #         # S <- S|set(set(A$vnames[i],A$vnames[j]))
#   #         
#   #         PA <- grp(A, disij, etype="directed", dir=-1, inclusive=FALSE)  # other parents
#   # 
#   #         for (k in seq_along(PA)){
#   #           # S <- S|set(set(A$vnames[i],A$vnames[j],A$vnames[PA[k]]))
#   #           S <- c(S, list(A$v[c(i,j,PA[k])])) #
#   #         }
#   #         for (k in seq_along(distail)){
#   #           # S <- S|set(set(A$vnames[i],A$vnames[j],A$vnames[distail[k]]))
#   #           S <- c(S, list(A$v[c(i,j,distail[k])])) #
#   #         }
#   #       }
#   #     }
#   #   }
#   # }
#   # #above "for" loop is to identify heads of two vertex and their tails
#   
#   # dists <- districts(A)
#   # lens <- lengths(dists)
#   # 
#   # for (d in which(lens >= 3)) {
#   #   di <- dists[[d]]
#   #   
#   #   for (i in di[seq_len(lens[d]-2)]) for (j in di[seq_len(lens[d]-i-1)+i]) {
#   #   
#   #   # di <- dis(A,A$v[i],sort = 2)
#   #   
#   #     if ((!(j %in% anc[[i]])) && (!(i %in% anc[[j]]))) {
#   #       for (k in di[seq_len(lens[d]-j)+j]) {
#   #   # for (j in seq_len(Lv-1-i)+i){
#   #     # if ((j %in% di) && (!(j %in% anc[[i]])) && (!(i %in% anc[[j]]))){
#   #       # for (k in seq_len(Lv-j)+j) {
#   #         if (!(k %in% c(anc[[i]], anc[[j]]))) {
#   #           B <- subGraph(A,c(anc[[i]], anc[[j]], anc[[k]]))
#   #           dai <- dis(B,i,sort = 2)
#   #           if ((j %in% dai) && (k %in% dai) &&
#   #               (!(i %in% anc[[k]])) && (!(j %in% anc[[k]]))){
#   #             # S <- S|set(set(A$vnames[i],A$vnames[j],A$vnames[k]))
#   #             S <- c(S, list(A$v[c(i,j,k)]))
#   #           }
#   #         }
#   #       }
#   #     }
#   #   }
#   # }
#   #   #Last "for" loop is to identify heads of three vertex
#   
#   val <- sapply(S, function(x) sum(2^x))
#   S <- S[!duplicated(val)]
#   val <- val[!duplicated(val)]
#   S <- S[order(val)]
#   S <- lapply(S, sort.default)
#   
#   return(S)
# }

##' Get the subset representation of an ordinary ADMG
##' 
##' @param graph ADMG
##' @param max_size largest set to consider
##' @param sort if 1, returns unique sets, if 2 returns sorted sets, if 3 returns sets ordered reverse lexicographically
##' @param r logical: use recursive heads? (Defaults to \code{FALSE})
##' 
##' @details This returns sets of the form \code{HuA}, where \code{A} is any 
##' subset of \code{tail(H)}.  This representation is shown to 
##' correspond precisely to ordinary Markov equivalence.
##' 
##' @export subsetRep
subsetRep <- function (graph, max_size, sort=1, r=FALSE) {
  if (missing(max_size)) max_size <- length(graph$v)
  
  out <- list()
  
  ## start with sets from undirected part of graph
  unG <- un(graph)
  und_comp <- graph[unG, etype="undirected"]
  if (any(lengths(sapply(graph$v, function (x) c(pa(graph, x), sib(graph, x)))) > 0 
          & (graph$v %in% unG))) stop("Not a euphonious/summary graph")
  
  clq <- cliques(und_comp)
  for (i in seq_along(clq)) {
    tmp <- powerSet(clq[[i]], m=max_size)[-1]
    out <- c(out, tmp)
  }

  ## remove redundant entries if required 
  if (sort > 0) {
    fn <- sapply(out, function(x) sum(2^(x-1)))
    if (any(duplicated(fn))) out <- out[!duplicated(fn)]
  }
  
  graph$edges$undirected <- NULL
  
  ## now obtain sets from directed part of graph
  ht <- headsTails(graph, r=r, max_head = max_size, by_district = FALSE)
  
  for (i in seq_along(ht$heads)) {
    if (all(ht$heads[[i]] %in% unG)) next
    ps <- powerSet(ht$tails[[i]], m=max_size-length(ht$heads[[i]]))
    tmp <- lapply(ps, function(x) c(ht$heads[[i]], x))
    out <- c(out, c(tmp))
  }
  
  if (sort > 1) {
    out <- lapply(out, sort.int)
  }
  if (sort > 2) {
    fn <- sapply(out, function(x) sum(2^(x-1)))
    if (any(duplicated(fn))) stop("There should not be a duplication, but something appears twice")
    out <- out[order(fn)]
  }
  
  out
}

he_val <- function (graph, skel=TRUE) {
  n <- nv(graph)
  if (nv(graph) > 7) stop("Only works for graphs with at most 7 vertices")
  ## also doesn't work if n=7 and skel=TRUE...
  ssr <- subsetRep(graph, max_size = 3, r = FALSE, sort = 2)
  if (skel) {
    ssr <- ssr[lengths(ssr) >= 2]
    return(sum(2^(setmatch(ssr, c(combn(n, 2, simplify = FALSE), 
                                  combn(n, 3, simplify = FALSE)))-1)))
  }
  else {
    ssr <- ssr[lengths(ssr) == 3]
    return(sum(2^(setmatch(ssr, combn(n, 3, simplify = FALSE))-1)))
  }
}

##' Get sets constrained by a conditional independence
##' 
##' @param ci a conditional independence or list of \code{ci} objects
##' @param check logical: should validity of \code{ci} be checked?
##' 
##' @export
constr_sets <- function (ci, sort=1, max=Inf, check=TRUE) {
  
  if (length(ci) == 0) return(integer(0))
  
  ## deal with a list of CIs first
  if (is.list(ci[[1]])) {
    return(lapply(ci, constr_sets, sort=sort, max=max, check=check))
  }
  
  ## make into a ci object
  if (class(ci) != "ci") ci <- as.ci(ci)
  ## check if requested
  if (check) {
    if (length(intersect(ci[[1]], ci[[2]])) > 0) stop("Independent sets should not overlap")
    ci[[1]] <- setdiff(ci[[1]], ci[[3]])
    ci[[2]] <- setdiff(ci[[2]], ci[[3]])
  }
  
  out <- list()
  if (is.finite(max)) sz <- list()
  ## use kronecker to obtain indices for subsets
  for (i in 1:3) {
    out[[i]] <- 0
    if (is.finite(max)) sz[[i]] <- 0
    
    for (j in seq_along(ci[[i]])) {
      ## add in representation of jth element of set
      out[[i]] <- kronecker(c(0,2^(ci[[i]][j] - 1)), out[[i]], "+")
      if (is.finite(max)) {
        sz[[i]] <- kronecker(0:1, sz[[i]], "+")
        out[[i]] <- out[[i]][sz[[i]] <= max]
        sz[[i]] <- sz[[i]][sz[[i]] <= max]
      }
    }
    if (i < 3) {
      out[[i]] <- out[[i]][-1]
      if (is.finite(max)) sz[[i]] <- sz[[i]][-1]
    }
  }
  
  ## combine to obtain representation of final sets
  out_full <- kronecker(out[[3]], kronecker(out[[2]], out[[1]], "+"), "+")
  if (is.finite(max)) {
    osz <- kronecker(sz[[3]], kronecker(sz[[2]], sz[[1]], "+"), "+")
    out_full <- out_full[osz <= max]
    osz <- osz[osz <= max]
  }
  
  if (sort > 0) {
    ## if requested, sort the output
    out_full <- sort.int(out_full)
  }
  
  return(out_full)
}
