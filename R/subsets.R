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
##' 
##' @details This returns sets of the form \code{HuA}, where \code{A} is any 
##' subset of \code{tail(H)}.  This representation is shown to 
##' correspond precisely to ordinary Markov equivalence.
##' 
##' @export subsetRep
subsetRep <- function (graph, max_size, sort=1) {
  if (missing(max_size)) max_size <- length(graph$v)
  
  ht <- MixedGraphs::headsTails(graph, r=FALSE, max_head = max_size, by_district = FALSE)
  
  out <- list()
  
  for (i in seq_along(ht$heads)) {
    ps <- powerSet(ht$tails[[i]], m=max_size-length(ht$heads[[i]]))
    tmp <- lapply(ps, function(x) c(ht$heads[[i]], x))
    out <- c(out, c(tmp))
  }
  
  if (sort > 1) {
    fn <- sapply(out, function(x) sum(2^(x-1)))
    if (any(duplicated(fn))) stop("There should not be a duplication, but something appears twice")
    out <- out[order(fn)]
    out <- lapply(out, sort.int)
  }
  
  out
}

