#' Generate a distribution from an ADMG model
#' 
#' Generate a distribution from an ADMG model
#' 
#' @param graph an object of class graph
#' @param dims integer vector of dimensions of distribution
#' @param map optionally, the output of \code{maps()}
#' @param r logical indicating whether or not recursive factorizations should
#' be used
#' @param alpha parameter to use in beta distribution
#' 
#' @return Object of class \code{moebius} giving generalized Moebius parameters
#' for a distribution in the model associated with \code{graph}.
#' 
#' @export
rADMGdist <- function(graph, dims, map, r=TRUE) {
  if (missing(graph)) graph <- rADMG(n=length(dims))
  if (missing(map)) map <- maps(graph, dims=dims, r=r)
  mobs <- mobs2 <- moebius(graph, dims=dims, r=r)
  ok <- FALSE
  gen <- function(x) {
    out <- rbeta(length(x), 2*x, 2*(1-x))
    dim(out) <- dim(x)
    out
  }
  while(!ok) {
    mobs2$q <- rapply(mobs$q, gen, how="replace")
    if (all(probdist(mobs2, map, graph) >= 0)) ok = TRUE
  }
  mobs2
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
