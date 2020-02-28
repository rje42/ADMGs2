#' Generate a distribution from an ADMG model
#' 
#' Generate a distribution from an ADMG model
#' 
#' @param graph an object of class graph
#' @param dims integer vector of dimensions of distribution
#' @param maps optionally, the output of \code{maps()}
#' @param r logical indicating whether or not recursive factorizations should
#' be used
#' @return Object of class \code{moebius} giving generalized Moebius parameters
#' for a distribution in the model associated with \code{graph}.
#' 
#' @export rADMG
rADMG <- function(graph, dims, maps, r=TRUE) {
  if (missing(maps)) maps <- maps(graph, dims=dims, r=r)
  mobs <- mobs2 <- moebius(graph, dims=dims, r=r)
  ok <- FALSE
  gen <- function(x) {
    out <- rbeta(length(x), 2*x, 2*(1-x))
    dim(out) <- dim(x)
    out
  }
  while(!ok) {
    mobs2$q <- rapply(mobs$q, gen, how="replace")
    if (all(probdist(mobs2, maps, graph) >= 0)) ok = TRUE
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
#' @param maps list containing maps for (use \code{maps}).
#' @param graph if maps not supplied, an object of class \code{graph}; use to
#' calculate maps.
#' @return A numeric array of dimensions \code{moebius$dims} giving the joint
#' distribution described by the Moebius parameters.
#' @author Robin Evans
#' @references Evans and Richardson.
#' 
probdist <-
function(moebius, maps, graph) {

  if (missing(maps)) maps = maps(graph, dims=moebius$dims, r=moebius$r)
  probs = rep(1, prod(maps$dim))
  q = moebius$q

  for (i in seq_along(maps$M)) {
    tmp = as.vector(maps$M[[i]] %*% exp(maps$P[[i]] %*% log(unlist(q[[i]]))))
    probs = probs*patternRepeat(tmp, maps$pa.dists[[i]], maps$dim)
  }

  array(probs, maps$dim)
}
