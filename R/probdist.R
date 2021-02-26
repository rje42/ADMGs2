#' Generate a distribution from an ADMG model
#' 
#' Generate a distribution from an ADMG model
#' 
#' @param graph an object of class graph
#' @param dims integer vector of dimensions of distribution
#' @param map optionally, the output of \code{maps()}
#' @param r logical indicating whether or not recursive factorizations should
#' be used
#' @param alpha parameter to scale Beta distribution down by. 
#' 
#' @details The random distribution is
#' obtained by starting with the uniform distribution, and adding an 
#' independent (scaled) normal random variables to each Moebius parameter.  We 
#' then scale the move to reach approximately the boundary of the simplex,
#' and then pick a value between the uniform distribution and the boundary 
#' point using a Beta distribution.  The parameter \code{alpha} being larger
#' favours distributions closer to the uniform distribution.
#' 
#' @return Object of class \code{moebius} giving generalized Moebius parameters
#' for a distribution in the model associated with \code{graph}.  
#' 
#' @export
rADMGdist <- function(graph, dims, map, r=TRUE, alpha=1) {
  if (missing(dims)) dims <- rep(2, nv(graph))
  if (missing(graph)) graph <- rADMG(n=length(dims))
  if (missing(map)) map <- maps(graph, dims=dims, r=r)
  mobs <- mobs2 <- moebius(graph, dims=dims, r=r)
  
  qs <- unlist(mobs$q)
  
  mv <- runif(length(qs), min=-3*qs/4, max=3*qs/4) + 
    rnorm(length(qs), sd=0.1/sqrt(qs))
  if (any(qs + mv < 0)) {
    message("Scaling down")
    mv <- mv/(-min(mv/qs))+1e-10
  }
  mobs2$q <- relist(qs+mv, mobs$q)
  
  fact_d <- 1
  fact_u <- 1
  
  ## if move unnecessarily small, then extend
  while (all(probdist(mobs2, map, graph) >= 0)) {
    fact_u <- fact_u*1.5
    if (any(qs + fact_u*mv < 0)) {
      fact_u <- fact_u/1.5
      break
    }
    mobs2$q <- relist(qs+mv*fact_u, mobs$q)
  }

  ## if move too large, then contract
  while (any(probdist(mobs2, map, graph) < 0)) {
    fact_d <- fact_d*1.5
    if (any(qs + mv/fact_d < 0)) {
      fact_d <- fact_d/1.5
      break
    }
    mobs2$q <- relist(qs+mv/fact_d, mobs$q)
  }
  
  ## now scale with a beta to account for dimension
  scal <- rbeta(1, length(qs)/alpha, 1)
  mobs2$q <- relist(qs+mv*fact_u/fact_d*scal, mobs$q)
  mobs2 <- c(list(p=probdist(mobs2, map, graph)), mobs2)
  if (!all(mobs2$p >= 0)) stop("Negative probability encountered")
  
  class(mobs2) <- "mparam"
  
  return(mobs2)
  # 
  # h1 <- which(lengths(mobs$heads) == 1)
  # 
  # ok <- FALSE
  # gen <- function(x) {
  #   len <- length(x)
  #   if (!is.null(dim(x))) full <- prod(dim(x) + 1)
  #   else full <- len + 1
  #   out <- c(rje::rdirichlet(1, rep(alpha, full)))
  #   out <- out[-(len+seq_len(full-len))]
  #   # out <- rbeta(length(x), 2*x*alpha, 2*(1-x)*alpha)
  #   dim(out) <- dim(x)
  #   out
  # }
  # mobs$q[h1] <- mobs2$q[h1] <- rapply(mobs$q[h1], gen, how="replace")
  # 
  # if (length(h1) < length(mobs$heads)) {
  #   while(!ok) {
  #     mobs2$q[-h1] <- rapply(mobs$q[-h1], gen, how="replace")
  #     if (all(probdist(mobs2, map, graph) >= 0)) ok = TRUE
  #   }
  # }
  # mobs2
}

#' Generate distributions from an ADMG model
#' 
#' Generate distributions from an ADMG model
#' 
#' @param n number of distributions to generate
#' @param graph an object of class graph
#' @param dims integer vector of dimensions of distribution
#' @param r logical indicating whether or not recursive factorizations should
#' be used
#' 
#' @details The random distribution is
#' obtained by starting with the uniform distribution, and adding an 
#' independent (scaled) normal random variables to each Moebius parameter.  We 
#' then scale the move to reach approximately the boundary of the simplex,
#' and then pick a value between the uniform distribution and the boundary 
#' point using a Beta distribution.  
#' 
#' @return Object of class \code{tables} giving distributions Markov with 
#' respect to a distribution in the model associated with \code{graph}.  
#' 
#' @importFrom contingency repTables
#' 
#' @export
rmixedgraph <- function(n, graph, dims, r=TRUE, alpha=1) {
  if (missing(dims)) dims <- rep(2, nv(graph))
  
  map <- maps(graph, dims=dims, r=r)

  ## produce tables object
  out <- repTables(n, function() rADMGdist(graph, dims, map, r=r, alpha=alpha)$p)
  
  return(out)
}


#' Calculate joint probabilities.
#' 
#' Given the generalized Moebius parametrization of an ADMG, calculates an
#' array of joint probabilities.
#' 
#' 
#' @param moebius Object of class \code{mparam}.
#' @param map list containing maps for (use \code{maps}).
#' @param graph if maps not supplied, an object of class \code{graph}; use to
#' calculate maps.
#' @return A numeric array of dimensions \code{moebius$dims} giving the joint
#' distribution described by the Moebius parameters.
#' @author Robin Evans
#' @references Evans and Richardson.
#' 
#' @export
probdist <-
function(moebius, map, graph) {

  if (missing(map)) map = maps(graph, dims=moebius$dim, r=moebius$r)
  probs = rep(1, prod(map$dim))
  q = moebius$q

  for (i in seq_along(map$M)) {
    tmp = as.vector(map$M[[i]] %*% exp(map$P[[i]] %*% log(unlist(q[[i]]))))
    probs = probs*patternRepeat(tmp, map$pa.dists[[i]], map$dim)
  }

  array(probs, map$dim)
}
