##' Function to generate random graphs
##' 
##' @param n number of vertices
##' @param N number of graph to simulate
##' @param p_dir,p_bi probability of a directed/bidirected edge
##' @param list logical: should output always be a list?
##' @param arid logical: should graph be made to be arid?
##' 
##' @details Uses a binomial number of edges based on \code{p_dir} and \code{p_bi}.
##' 
##' @export
rADMG <- function(n, N=1, p_dir=1/n, p_bi=1/n, list=FALSE, arid=FALSE) {
  ne <- choose(n, 2)

  # requireNamespace(Matrix)
  
  if (N > 1 || list) out2 <- list()
  
  for (i in seq_len(N)) {
    ## sample edges to include  
    nedg <- rbinom(2, size=ne, prob=c(p_dir, p_bi))

    # dir <- getEdges(n=n, nedge=nedg[1])
    dir_ind <- sample(ne, nedg[1], replace=FALSE)
    bi_ind <- sample(ne, nedg[2], replace=FALSE)

    nv <- c(0,cumsum(seq_len(n-1)))
        
    dir <- bi <- matrix(nrow=2, ncol=0)
    # done_dir <- rep(FALSE, nedg[1])
    # done_bi <- rep(FALSE, nedg[2])
    # new <- rep(TRUE, nedg[1])
    # 
    # for (i in seq_len(n)[-1]) {
    #   done_dir <- dir_ind < i
    #   if (any(new & done_dir)) {
    #     dir <- cbind(dir, rbind(dir_ind[new & done_dir], rep(i, sum(new & done_dir))))
    #     new[done_dir] <- FALSE
    #   }
    #   dir_ind[new] <- dir_ind[new] - i
    # }
    
    # for (i in seq_len(n)[-1]) {
    #   new <- which(!done_dir & dir_ind <= nv[i])
    #   dir <- cbind(dir, matrix(c(dir_ind[new] - nv[i-1],rep(i, length(new))), nrow=2, byrow=TRUE))
    #   done_dir[new] = TRUE
    # 
    #   new <- which(!done_bi & dir_ind <= nv[i])
    #   bi <- cbind(bi, matrix(c(bi_ind[new] - nv[i-1],rep(i, length(new))), nrow=2, byrow=TRUE))
    #   done_bi[new] = TRUE
    # 
    #   # printCount(i, n, )
    # }

    # class(dir) <- class(bi) <- "edgeMatrix"
    
    
    nv <- cumsum(seq_len(n-1))
            
    if (p_dir > 0) {
      ## collect directed edges
      out <- .rowSums(outer(dir_ind, nv, `>`), nedg[1], n-1)+2
      dir <- matrix(c(dir_ind - c(0,nv)[out-1], out), nrow=2, byrow=TRUE)
      names(dir) <- NULL
      class(dir) <- "edgeMatrix"
    }

    if (p_bi > 0) {
      ## collect bidirected edges
      out <- rowSums(outer(bi_ind, nv, `>`))+2
      bi <- rbind(bi_ind - c(0,nv)[out-1], out)
      names(bi) <- NULL
      class(bi) <- "edgeMatrix"
    }
    
    ## put together edges
    if (p_bi > 0 && p_dir > 0) edg <- list(directed=dir, bidirected=bi)
    else if (p_bi > 0) edg <- list(bidirected=bi)
    else if (p_dir > 0) edg <- list(directed=dir)
    else edg <- list()
    class(edg) <- "edgeList"
    
    out <- mixedgraph(n, edges=edg)
    
    if (arid) {
      out <- aridProj(out, maximal=TRUE)
    }
    if (N > 1 || list) out2[[i]] <- out
    else out2 <- out
  }
  
  return(out2)
}

# 
# getEdges <- function(n, nedge) {
#   nv <- c(0,cumsum(seq_len(n-1)))
#   
#   
# }