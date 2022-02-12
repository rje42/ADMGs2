##' Turn list into a conditional independence
##' 
##' @param x either a list with two or three elements, or a vector
##' @param ... if \code{x} not a list, one or two more arguments
##' 
##' @examples 
##' as.ci(list(1,2))
##' as.ci(list(1,2,3))
##' as.ci(1,2)
##' as.ci(1,2,3)
##' 
##' @export
as.ci <- function(x, ...) {
  if (!is.list(x)) {
    args <- list(...)
    if (length(args) == 1) x <- list(x, args[[1]])
    else if (length(args) == 2) x <- list(x, args[[1]], args[[2]])
    else stop("must provide 2 or 3 arguments or first argument must be a list")
  }

  if (length(x) == 2) {
    x[[3]] <- integer(0)
  }
  else if (!(length(x) == 3)) stop("should have length 2 or 3")
  
  x <- lapply(x, as.integer)
  class(x) <- "ci"
  
  x
}

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
  if (!is_ADMG(graph)) stop("Input should be an acyclic directed mixed graph")
  
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
    out <- split_ci(out)
  }
  
  standardize_cis(out)
}

##' Operations for conditional independences
##' 
##' Operations on conditional independence objects
##' 
##' @param cis object of class \code{ci} or list thereof
##' 
##' @return A list with the relevant output in, one entry for each conditional 
##' independence put into \code{cis}. 
##' 
##' @name ci_operations
NULL

##' @describeIn ci_operations Split conditional independences into univariate pieces
##' @export
split_ci <- function (cis) {
  if ("ci" %in% class(cis)) cis <- list(cis)
  
  new_indeps <- list()
  rmv <- integer(0)
  
  ## start by splitting any second terms
  for (i in seq_along(cis)) {
    if (length(cis[[i]][[2]]) > 1) {
      vnew_indeps <- vector(mode="list", length=length(cis[[i]][[2]]))
      cond <- cis[[i]][[3]]
      for (j in seq_along(cis[[i]][[2]])) {
        vnew_indeps[[j]] <- list(cis[[i]][[1]], cis[[i]][[2]][j], cond)
        cond <- c(cond, cis[[i]][[2]][j])
        class(vnew_indeps[[j]]) <- "ci"
      }
      rmv <- c(rmv, i)
      new_indeps <- c(new_indeps, vnew_indeps)
    }
  }
  if (length(rmv) > 0) cis <- c(cis[-rmv], new_indeps)

  ## reset new_indeps and rmv
  new_indeps <- list()
  rmv <- integer(0)
  
  ## now split any first terms  
  for (i in seq_along(cis)) {
    if (length(cis[[i]][[1]]) > 1) {
      vnew_indeps <- vector(mode="list", length=length(cis[[i]][[1]]))
      cond <- cis[[i]][[3]]
      for (j in seq_along(cis[[i]][[1]])) {
        vnew_indeps[[j]] <- list(cis[[i]][[1]][j], cis[[i]][[2]], cond)
        cond <- c(cond, cis[[i]][[1]][j])
        class(vnew_indeps[[j]]) <- "ci"
      }
      rmv <- c(rmv, i)
      new_indeps <- c(new_indeps, vnew_indeps)
    }
  }
  if (length(rmv) > 0) cis <- c(cis[-rmv], new_indeps)
  
  
  cis <- rapply(cis, as.integer, how="replace")
  cis <- lapply(cis, function(x) `class<-`(x, "ci"))

  cis
}

##' @describeIn ci_operations Standardize conditional independences
##' @export
standardize_cis <- function (cis) {
  if ("ci" %in% class(cis)) cis <- list(cis)
  
  cis2 <- purrr::transpose(cis)
  
  cis2[[1]] <- mapply(function(x,y) sort.int(unique.default(setdiff(x,y))), 
                      cis2[[1]], cis2[[3]], SIMPLIFY = FALSE)
  cis2[[2]] <- mapply(function(x,y) sort.int(unique.default(setdiff(x,y))), 
                      cis2[[2]], cis2[[3]], SIMPLIFY = FALSE)
  cis2[[3]] <- lapply(cis2[[3]], function(x) sort.int(unique.default(x)))
  
  rmv <- (lengths(cis2[[1]])==0) | (lengths(cis2[[2]])==0)
  cis2 <- purrr::transpose(purrr::transpose(cis2)[!rmv])
  if (length(cis2) == 0) return(list())
  
  if (any(lengths(mapply(intersect, cis2[[1]], cis2[[2]])) > 0)) stop("Should not be overlap between two main sets")
  
  minA <- sapply(cis2[[1]], min)
  minB <- sapply(cis2[[2]], min)
  tmpA <- tmpB <- vector("list", length(cis2[[1]]))
  tmpA[minA < minB] <- cis2[[2]][minA < minB]
  tmpA[minA > minB] <- cis2[[1]][minA > minB]
  tmpB[minA < minB] <- cis2[[1]][minA < minB]
  tmpB[minA > minB] <- cis2[[2]][minA > minB]
  
  cis2[[1]] <- tmpA
  cis2[[2]] <- tmpB
  
  cis <- purrr::transpose(cis2)
  cis <- rapply(cis, as.integer, how="replace")
  cis <- lapply(cis, function(x) `class<-`(x, "ci"))
  

  cis
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
