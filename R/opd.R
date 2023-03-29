##' Get the ordered power DAG
##'
##' Returns the ordered power DAG for an ordinary Markov ADMG
##'
##' @param graph an ADMG of class \code{mixedgraph}
##' @param topOrd a topological ordering of the vertices in \code{graph}
## @param r logical: recursive formulation?
##' @param v optionally the vertices to obtain components for
##' @param include_empty logical: should independences be included even if there are no variables in one set?
##' 
##' @details Returns the ordered power DAG, as defined in Hu and Evans (2022),
##' giving independences implied by the local Markov property and sufficient to
##' define the model.
##' 
##' @references Z. Hu and R.J. Evans. Towards standard imsets for maximal ancestral
##' graphs, _arXiv:2208.10436_, 2022.
##' 
##' @export
opd <- function (graph, topOrd, v, include_empty=FALSE) {
  if (missing(topOrd)) topOrd <- topologicalOrder(graph)
  else {
    if (!is.subset(graph$v, topOrd)) stop("Not a valid topological ordering")
    topOrd <- topOrd[!duplicated(topOrd)]
  }
  if (missing(v)) v <- graph$v
  
  # graph2 <- graph[topOrd, order=TRUE]
  graph2 <- graph
  PA <- lapply(seq_along(graph$vnames), function(x) pa(graph, x))
  n <- length(topOrd)
  
  # v2 <- match(v, topOrd)
  out <- mixedgraph(n=0)
  indep <- list()
  
  for (i in seq_len(n)) {
    vi <- topOrd[i]
    ## get maximal intrinsic set for vertex i
    mdis <- dis(graph2[topOrd[seq_len(i)]], vi)
    padis <- unique.default(c(mdis, unlist(PA[mdis])))
    gr2 <- graph2[padis]

    # ## add code if recursive part implemented
    # if (r) {
    #   ## code to sever edges into parents variables
    # }
    
    ## add global independence
    if (length(padis) < i || include_empty) {
      indep <- c(indep, list(as.ci(list(vi, setdiff(topOrd[seq_len(i)], padis), setdiff(padis, vi)))))
    }
    
    ## add to list
    curr <- list(gr2)
    curr2 <- list()
    sec <- list(barren(gr2))
    ct <- 0
    pa <- NA
    
    new <- TRUE
    
    while (length(curr) > 0) {
      for (k in seq_along(curr)) {
        ## now go through each new graph and remove each barren vertex (other than vi)
        for (l in setdiff(rev(barren(curr[[k]])), vi)) {
          poss <- dis(curr[[k]][-l], vi, sort=2)
          ## if we already have it, then skip
          if (setmatch(list(poss), sec, nomatch = 0L) > 0) next

          ## otherwise add it to the list, and record the parent set
          sec <- c(sec, list(barren(curr[[k]], poss)))
          poss_all <- c(poss, setdiff(pa(curr[[k]], poss), poss))
          curr2 <- c(curr2, list(curr[[k]][poss_all]))
          pa <- c(pa, ct + k)
          
          v_lost <- setdiff(setdiff(curr[[k]]$v, l), last(curr2)[[1]]$v)
          if (include_empty || length(v_lost) > 0) {
            poss_all2 <- setdiff(poss_all, vi)
            indep <- c(indep, list(as.ci(vi, v_lost, setdiff(poss_all2, v_lost))))
          }
        }
      }
      ## having finished this level, set the curr to curr2 and reset curr2 to be empty
      ct <- ct + length(curr)
      curr <- curr2
      curr2 <- list()
    }
    
    if(n > 9) vnms <- sapply(sec, function(x) paste0("s",paste(rev(x), collapse=",")))
    else vnms <- sapply(sec, function(x) paste0("s",paste(rev(x), collapse="")))
    # browser()
    nvs <- nv(out)
    out <- MixedGraphs::addNodes(out, length(sec), vnames = vnms)
    if (length(sec) > 1) {
      dedges <- apply(nvs + cbind(pa, seq_along(sec))[-1,,drop=FALSE], 1, c, simplify = FALSE)
      out <- addEdges(out, makeEdgeList(directed=eList(dedges)))
    }
  }
  
  attr(out, "indep") <- indep
  out
}

##' Obtain the refined ordered Markov property
##' 
##' @param graph an ADMG of class \code{mixedgraph}
##' @param topOrd optional topological ordering
##' @param power_DAG a power DAG for \code{graph}
##' 
##' @details Computes a power DAG and extracts the \code{indep} attribute.
##' 
##' @export
romp <- function (graph, topOrd, power_DAG, group=FALSE) {
  if (missing(power_DAG)) power_DAG <- opd(graph, topOrd)
  
  out <- attr(power_DAG, "indep")
  if (group) out <- group_indeps(out)
  
  return(out)
}

