#' Fit data to an ADMG model
#' 
#' Fit discrete data to the Markov structure implied by an acyclic directed
#' mixed graph.
#' 
#' Fits data using coordinate-wise block descent algorithm, with an Armijo line
#' search.  Details in Evans and Richardson (2010).
#' 
#' Recursive factorizations allow for modelling of Verma-constraints.  See
#' Richardson et al (2012).
#' 
#' @aliases fitADMG print.mixed_fit
#' @param dat The data, as an array of counts or a data frame whose final
#' column contains the counts.
#' @param graph An ADMG, as an object of class \code{graph}.
#' @param r Logical indicating whether or not recursive factorizations should
#' be used.
#' @param tol Numeric; if log-likelihood increases by less that \code{tol} in
#' one step, procedure stops.
#' @param SEs logical: should standard errors be calculated?
#' @param sparse Should sparse matrices be used?  Requires package
#' \code{Matrix}.
#' @param quietly Logical indicating whether output should be suppressed.
#' @param use_optim should \code{optim} be used for fitting?
#' 
#' @return An object of class \code{mixed_fit}.  This is a list containing
#' (amongst other things): \item{params}{An object of class \code{Mparams},
#' containing the Moebius parameters which maximise the likelihood.}
#' \item{ll}{Value of the log-likelihood at the maximum.}
#' @note %% ~~further notes~~
#' @section Warning: For the algorithm to be guaranteed to work correctly, all
#' counts for marginal tables consisting of districts and their parents should
#' be positive.  A warning will be produced if this is not so.
#' @author Robin Evans, Ilya Shpitser
#' @seealso \code{\link{summary.mixed_fit}} \code{\link{autoFit}}
#' @references Evans, R.J. and Richardson, T.S. (2010) - Fitting acyclic
#' directed mixed graphs to binary data. \emph{UAI-10}.
#' 
#' Richardson, T.S., Robins, J.M., Shpitser, I., Evans, R.J. (2012) -
#' @keywords graphs
#' @examples
#' 
#' data(gss_small)
#' data(gr1, package="MixedGraphs")
#' 
#' out = fitADMG(gss_small, gr1)
#' summary(out)
#' # not a good fit!
#' 
#' @export fitADMG
fitADMG <- function(dat, graph, r = TRUE, tol = sqrt(.Machine$double.eps), 
                    SEs = TRUE, sparse = FALSE, quietly = TRUE, use_optim=TRUE) {
  n = max(graph$v)
  if (tol < 0) stop("Tolerance must be positive")
  
  if (is.array(dat)) dims = dim(dat)
  else dims = sapply(dat[,-ncol(dat),drop=FALSE], max)+1
  
  params = moebius(graph, dims=dims, r=r)
  if (!quietly) cat("Getting maps\n")
  maps0 = maps(graph, sparse = sparse, dims = dims, r=r)
  
  ## # ORDER DATA AND ADD 0s IF NECESSARY
  ## dat = orderData(dat)
  ## new.dat = cbind(combinations(dims), 0)
  ## rows = prodlim::row.match(dat[,seq_len(n)], new.dat[,seq_len(n)])
  ## new.dat[rows, n+1] = dat[, n+1]
  ## dat = array(new.dat[,n+1], dim=dims)
  
  # SUFFICIENT STATISTICS ARE COUNTS OVER D \cup pa(D), FOR DISTRICT D.
  # WRITE AS A LIST OF THESE.
  # if(is.data.frame(dat)) {
  #   dat = plyr::daply(dat, seq_len(n), function(x) x[[n+1]])
  # }
  
  dist.dat = list()
  for (d in seq_along(maps0$dists)) {
    if (is.array(dat)) dist.dat[[d]] = margin.table(dat, maps0$pa.dists[[d]])
    else {
      comb <- rje::combinations(dims[maps0$pa.dists[[d]]])
      dist.dat[[d]] <- array(dims[maps0$pa.dists[[d]]])

      for (i in seq_len(nrow(comb))) {
        inc <- (rowSums(dat[,maps0$pa.dists[[d]],drop=FALSE] == comb[rep(i, nrow(dat)),,drop=FALSE]) == length(maps0$pa.dists[[d]]))
        dist.dat[[d]][i] <- sum(dat[inc, n+1])
      }
    }
    if (isTRUE(any(dist.dat[[d]] == 0))) warning("Zero counts in sufficient statistics")
  }

  # INITIAL TOLERANCE FOR .doone
  tol2 = 1e-1
  move = 2*tol

  if (!quietly) cat("Fitting: ")

  if (use_optim) {
    out_optim <- suppressWarnings(optim(unlist(params$q), fn = function(x) -.logLik(dat=dist.dat, q=relist(x, params$q), maps=maps0),
                       control=list(reltol=100*tol)))
    params$q <- relist(out_optim$par, params$q)

    # now do Newton's method to get even closer
    out_optim <- suppressWarnings(optim(unlist(params$q), fn = function(x) -.logLik(dat=dist.dat, q=relist(x, params$q), maps=maps0),
                       gr = function(x) -.DlogLik(dat=dist.dat, q=relist(x, params$q), maps=maps0),
                       method="BFGS", control=list(reltol=tol)))
    params$q <- relist(out_optim$par, params$q)
    
    # print(out_optim$value)
    # print((function(x) -.logLik(dat=dist.dat, q=relist(x, params$q), maps=maps0))(out_optim$par))
    # print(numDeriv::grad(out_optim$par, func=function(x) -.logLik(dat=dist.dat, q=relist(x, params$q), maps=maps0)))
  }
  else {
    count = 1
    
    while (move > tol || tol2 > sqrt(.Machine$double.eps)) {
      new.params = .doone(dist.dat, graph, params, maps0, tol = tol2, quietly = quietly)
      move = .logLik(dist.dat, new.params$q, maps0) - .logLik(dist.dat, params$q, maps0)
      params = new.params
      tol2 = tol2/10
      if (!quietly) cat(count," ")
      # print(params$q[[1]][[15]][[1]][13:16])
      count = count+1
    }
    if (!quietly) cat("\n")
  }
  
  out = list(params = params, ll = .logLik(dist.dat, params$q, maps0), dat=dat, graph=graph, dim=dims, maps=maps0, r=r)
  class(out) = "mixed_fit"
  
  if (SEs) {
    requireNamespace("numDeriv")
    if (!quietly) cat("Computing SEs")
    
    if (prod(dims) > 5e3) {
      f <- function(x) {
        q <- relist(x, params$q)
        -.logLik(dist.dat, q, maps0)
      }
      FIM <- numDeriv::hessian(f, x=unlist(params$q))
    }
    else FIM <- .fisher_mixed_fit(out)

    out$SEs <- sqrt(diag(solve.default(FIM)))
        
    if (!quietly) cat("\n")  
  }
  else out$SEs <- NULL

  out
}

