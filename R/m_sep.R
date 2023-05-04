##' Find separator sets
##' 
##' Implement the algorithms of van der Zander et al. (2014) to find separator sets.
##' 
##' @param graph a summary graph of class \code{mixedgraph}
##' @param X,Y sets to separate
##' @param Z set to include in conditioning set
##' @param R maximal set to look in
##' @param check logical: check that the sets are suitably disjoint?
##' 
##' @details \code{findSep} finds a set that m-separates \code{X} and \code{Y}, 
##' that must contain vertices in \code{Z} and may contain vertices in \code{R}.  
##' \code{listSep} exhaustively lists all such sets.  
##' 
##' \code{allSep} returns the maximal set m-separated from \code{X} by 
##' \code{Z}. \code{whSep} returns the subset of \code{Z} that is m-separated 
##' from \code{X} by the remainder of \code{Z}.
##' 
##' The algorithms for \code{findSep} and \code{listSep} are taken
##' from van der Zander et al. (2014).  
##' 
##' @references van der Zander et al. Constructing Separators and Adjustment Sets in Ancestral Graphs.
##' \emph{UAI Workshop on Causal Inference: Learning and Prediction}, 2014.
##' 
##' @export
findSep <- function (graph, X, Y, Z=integer(0), R, check=TRUE) {
  ## initialize variables
  if (missing(R)) {
    R <- setdiff(graph$v, c(X,Y))
  }
  if (check) {
    if (length(intersect(X,Y)) > 0) return(NA)
    X <- setdiff(X, Z)
    Y <- setdiff(Y, Z)
    Z <- intersect(Z, R)
    R <- union(setdiff(R, c(X,Y)), Z)
  }
  
  ## take the intersection of the anterior over X,Y,Z
  prop_C <- intersect(ant(graph, c(X,Y,Z)), setdiff(R, Z))
  ## return this set if and only if the m-separation holds
  if (m_sep(graph, X, Y, c(Z,prop_C))) return(c(Z,prop_C))
  else return(NA)
}

##' @describeIn findSep list all m-separations
##' @export
listSep <- function (graph, X, Y, Z=integer(0), R, check=TRUE) {
  ## initialize variables
  if (missing(R)) {
    R <- setdiff(graph$v, c(X,Y))
  }
  if (check) {
    if (length(intersect(X,Y)) > 0) return(NA)
    X <- setdiff(X, Z)
    Y <- setdiff(Y, Z)
    Z <- intersect(Z, R)
    R <- union(setdiff(R, c(X,Y)), Z)
  }

  ## check the separation with the given Z and R  
  out <- list(findSep(graph, X, Y, Z, R, check=FALSE))
  if (!isTRUE(is.na(out[[1]]))) {
    ## if there is one, then check for more
    RZd <- setdiff(R, Z)

    if (length(RZd) > 0) {
      ## if there is a gap between Z and R, then recall the function
      out <- c(listSep(graph, X, Y, Z, setdiff(R, RZd[1]), check=FALSE),
               listSep(graph, X, Y, c(Z, RZd[1]), R, check=FALSE))
    }
  }
  else out <- list()
  
  out
}


##' @describeIn findSep give maximal m-separated set
##' @export
allSep <- function (graph, X, Z=integer(0), R, check=TRUE) {
  ## initialize variables
  if (missing(R)) {
    R <- setdiff(graph$v, X)
  }
  if (check) {
    X <- setdiff(X, Z)
    # Y <- setdiff(Y, Z)
    Z <- intersect(Z, R)
    R <- union(setdiff(R, X), Z)
  }
  
  ## get set of things that are independent of X given Z
  gr_m <- moralize(graph, c(X,Z))
  C <- setdiff(gr_m$v, grp(gr_m[-Z], X))
  C <- setdiff(C, pathConnected(gr_m, X, Z, dir=0)) # remove not m-separated things in Z
  
  return(intersect(C, R))
}

##' @describeIn findSep find m-separated subset
##' @export
whSep <- function (graph, X, Z, check=TRUE) {
  ## initialize variables
  if (check) {
    X <- setdiff(X, Z)
  }
  
  gr_m <- moralize(graph, c(X,Z))
  out <- setdiff(Z, pathConnected(gr_m[ant(gr_m, c(X,Z))], X, Z, dir=0))  ## maybe don't need the ant() here...
  
  return(out)
}
