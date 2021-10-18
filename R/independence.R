##' Define conditional independence class
##' 
##' @exportS3Method print ci
print.ci <- function(x, ...) {
  cat(x[[1]], sep=", ")
  cat(" _||_ ")
  cat(x[[2]], sep=", ")
  if (length(x[[3]]) > 0) {
    cat(" | ")
    cat(x[[3]], sep=", ")
  }
  cat("\n")
  invisible(x)
}

##' Local Markov property for ADMGs
##' 
##' Get a list of independences implied by the local Markov property applied to
##' an ADMG.
##' 
##' @param graph an ADMG as an object of class \code{mixedgraph}
##' @param unique logical: should overlapping independences be removed?
##' @param split logical: should independences involve only singletons?
##' 
##' @export
localMarkovProperty <- function (graph, unique=TRUE, split=FALSE) {
  if (!is.mixedgraph(graph)) stop("Input must be a mixed graph object")
  if (!is.ADMG(graph)) stop("Input should be an acyclic directed mixed graph")
  
  ancSet <- anSets(graph, sort=3)
  ord <- topologicalOrder(graph)
  
  out <- list()
  
  ## get list of all independences
  for (i in seq_along(ancSet)) {
    idx <- which.max(match(ancSet[[i]], ord))
    v <- ancSet[[i]][idx]
    mbl <- setdiff(mb(graph[ancSet[[i]]], v, sort=2), v)
    if (length(c(v,mbl)) < length(ancSet[[i]])) {
      resid <- setdiff(ancSet[[i]], c(mbl, v))
      ci <- list(v,resid, mbl)
      class(ci) <- "ci"
      
      ## if required, remove anything implied by new independence
      if (unique) {
        rmv <- integer(0)
        for (j in seq_along(out)) {
          if (out[[j]][[1]] == v && indep_implies(ci, out[[j]])) rmv <- c(rmv, j)
        }
        if (length(rmv) > 0) out <- out[-rmv]
      }
      
      out <- c(out, list(ci))
    }
  }
  
  if (unique) {
    out <- out[!duplicated(out)]
  }
  
  ## if required, split into univariate independences
  if (split) {
    new_indeps <- list()
    rmv <- integer(0)
    
    for (i in seq_along(out)) {
      if (length(out[[i]][[2]]) > 1) {
        vnew_indeps <- vector(mode="list", length=length(out[[i]][[2]]))
        cond <- out[[i]][[3]]
        for (j in seq_along(out[[i]][[2]])) {
          vnew_indeps[[j]] <- list(out[[i]][[1]], out[[i]][[2]][j], cond)
          cond <- c(cond, out[[i]][[2]][j])
          class(vnew_indeps[[j]]) <- "ci"
        }
        rmv <- c(rmv, i)
        new_indeps <- c(new_indeps, vnew_indeps)
      }
    }
    out <- c(out[-rmv], new_indeps)
  }
  
  out
}

## Tests if one independence implies another via graphoids 2 and 3
indep_implies <- function(I1, I2) {
  ul1 <- unlist(I1)
  ul2 <- unlist(I2)
  
  ## swap around sets if necessary
  if (!(is.subset(I2[[1]], I1[[1]]))) I2 <- I2[c(2,1,3)]

  ## variables in I2 must be present in I1  
  if (!is.subset(ul2, ul1)) return(FALSE)
  if (!(is.subset(I2[[1]], I1[[1]]) && is.subset(I2[[2]], I1[[2]])) &&
      !(is.subset(I2[[1]], I1[[1]]) && is.subset(I2[[2]], I1[[2]]))) return(FALSE)
  
  ## then just check that all variables are present somewhere in I2
  if (!is.subset(I1[[3]], ul2)) return(FALSE)
  
  return(TRUE)
}

# ## Remove redundancy in list of independences
# rmvRed <- function(I) {
#   ulI <- lapply(I, function(x) sort.int(unlist(x)))
#   
#   
#   # 
#   # ## swap around sets if necessary
#   # if (!(is.subset(I2[[1]], I1[[1]]))) I2 <- I2[c(2,1,3)]
#   # 
#   # ## variables in I2 must be present in I1  
#   # if (!is.subset(ul2, ul1)) return(FALSE)
#   # if (!(is.subset(I2[[1]], I1[[1]]) && is.subset(I2[[2]], I1[[2]])) &&
#   #     !(is.subset(I2[[1]], I1[[1]]) && is.subset(I2[[2]], I1[[2]]))) return(FALSE)
#   
#   return(TRUE)
# }
