#' Fit data to an ADMG model
#' 
#' Fit discrete data to the Markov structure implied by an acyclic directed
#' mixed graph or summary graph.
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
#' @param graph An ADMG (or summary graph), as an object of class \code{graph}.
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
#' 
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
#' @export
fitADMG <- function(dat, graph, r = TRUE, tol = sqrt(.Machine$double.eps), 
                    SEs = TRUE, sparse = FALSE, quietly = TRUE, use_optim=TRUE) {
  n = length(graph$vnames)
  if (tol < 0) stop("Tolerance must be positive")
  
  if (is.array(dat)) dims = dim(dat)
  else if (is.data.frame(dat)) dims = sapply(dat[,-ncol(dat),drop=FALSE], max)+1
  else stop("invalid format for data: should be an array or data.frame object")
  
  ## remove any undirected part
  U <- un(graph)
  
  if (length(U) > 0) {
    gr_un <- graph[U, etype="undirected"]
    
    if (is.array(dat)) dat_u <- margin.table(dat, U)
    else {
      comb <- unique.data.frame(dat[,U])
      
      dat_u <- cbind(comb, count=NA)
      
      for (i in seq_len(nrow(comb))) {
        inc <- (rowSums(dat[,U,drop=FALSE] == comb[rep(i, nrow(dat)),,drop=FALSE]) == length(U))
        dat_u[i,ncol(dat_u)] <- sum(dat[inc, n+1])
      }
    }
    
    out_ug <- fitUG(dat_u, gr_un)
    
    if (length(U) == nv(graph)) return(out_ug)
  }
  else out_ug <- list(ll=0, p=0, df=0)
  
  graph <- graph[etype=c("directed", "bidirected")]
  graph <- fix(graph, U)
  
  ## get maps for ADMG component
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
  
  ## get sufficient statistics over districts
  dist.dat = list()
  for (d in seq_along(maps0$dists)) {
    if (is.array(dat)) dist.dat[[d]] = margin.table(dat, maps0$pa.dists[[d]])
    else if (is.data.frame(dat)) {
      comb <- rje::combinations(dims[maps0$pa.dists[[d]]])
      dist.dat[[d]] <- array(dims[maps0$pa.dists[[d]]])

      ## check the count for each district
      for (i in seq_len(nrow(comb))) {
        inc <- (rowSums(dat[,maps0$pa.dists[[d]],drop=FALSE] == comb[rep(i, nrow(dat)),,drop=FALSE]) == 
                  length(maps0$pa.dists[[d]]))
        dist.dat[[d]][i] <- sum(dat[inc, n+1])
      }
    }
    else stop("invalid format for data: should be an array or data.frame object")
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
  
  out = list(params = params, ll = .logLik(dist.dat, params$q, maps0) + out_ug$ll, 
             p = length(unlist(params$q)) + out_ug$p, dat=dat, graph=graph, 
             dim=dims, maps=maps0, r=r)
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

#' Fit data to an undirected graphical model
#' 
#' Fit discrete data to the Markov structure implied by an undirected graph.
#' 
#' @param dat The data, as an array of counts or a data frame whose final
#' column contains the counts.
#' @param graph An ADMG, as an object of class \code{graph}.
#' @param tol Numeric; if log-likelihood increases by less that \code{tol} in
#' one step, procedure stops.
#' @param SEs logical: should standard errors be calculated?
#' @param quietly Logical indicating whether output should be suppressed.
#' 
#' @return An object of class \code{u_fit}.  This is a list containing
#' (amongst other things): 
#' \item{ll}{Value of the log-likelihood at the maximum.}
#' 
#' @section Warning: For the algorithm to be guaranteed to work correctly, all
#' counts for marginal tables consisting of cliques and their parents should
#' be positive.  A warning will be produced if this is not so.
#' @author Robin Evans
#' @references Lauritzen (1996), Graphical Models, OUP.
#' @keywords graphs
#' 
#' @import dmod
#' @export
fitUG <- function(dat, graph, tol = sqrt(.Machine$double.eps), SEs = TRUE, 
                  quietly = TRUE) {
  if (!is.UG(graph)) stop("Must be an undirected graph")
  
  # A <- convert(graph, "graphNEL")

  ## get data into array format
  if (is.array(dat)) {
    dims <- dim(dat)
    dat2 <- as.table(dat)
  }
  else if (is.data.frame(dat)) {
    n <- ncol(dat) - 1
    dims <- sapply(dat[-(n+1)], max) + 1
    dat2 <- array(0, dim=dims)
    key <- c(1L, cumprod(dims)[-n])
    
    for (i in seq_len(nrow(dat))) {
      wh <- sum(dat[i,-(n+1)]*key) + 1
      dat2[wh] <- dat2[wh] + dat[i,n+1]
    }
    dimnames(dat2) <- lapply(dims, function(x) seq_len(x) - 1)
    # names(dimnames(dat2)) <- colnames(dat)[-(n+1)]

    dat2 <- as.table(dat2)
  }
  else stop("invalid format for data: should be an array or data.frame object")
  
  lls <- sum(dat2[dat2 > 0]*log(dat2[dat2 > 0]/sum(dat2)))
  
  if (length(dims) != nv(graph)) stop("Graph size does not match dimension of data")
  
  ### get formula
  clq <- cliques(graph)
  # form <- paste0("~ ", paste(sapply(clq, function(x) 
  #   paste(names(dat)[x], collapse=":")), collapse=" + "))
  # form <- as.formula(form)
  
  out0 <- loglin(dat2, clq, fit=TRUE, eps=sqrt(.Machine$double.eps), print=FALSE)
  out <- list(dev=out0$lrt, ll=lls-out0$lrt/2, lls=lls, fit=out0$fit, 
              p=prod(dims)-1-out0$df, df=out0$df)
  
  # out2 <- gRim::dmod(formula=form, data=dat2)
  
  # ## get saturated log-likelihood
  # cts_nz <- cts[cts > 0]
  # n_dat <- sum(cts_nz)
  # lls <- sum(cts_nz*log(cts_nz/n_dat))
  # 
  # ## get cliques, and consequent sufficient statistics
  # clqs <- cliques(graph, sort = 2)
  # tab <- vector(mode="list", length=length(clqs))
  # 
  # for (i in seq_along(clqs)) {
  #   if (is.array(dat)) tab[[i]] <- margin.table(dat, clqs[[i]])
  #   else {
  #     comb <- rje::combinations(dims[clqs[[i]]])
  #     for (j in seq_along(nrow(comb))) {
  #       wh <- (dat[,clqs[[i]]] == rep(comb[j,], each=nrow(dat)))
  #       tab[[i]][j] <- sum(dat[wh,ncol(dat)])
  #     }
  #   }
  # 
  #   
  # }  
  
  out
}


# ## run IPF
# run_ipf <- function(tab, clqs, control=list(), verbose=FALSE) {
# 
#   if (length(tab) != length(clqs)) stop("Must have same number of cliques as tables")
#     
#   # GET CONTROL PARAMETERS OR USE DEFAULTS
#   con = list(ep = .Machine$double.eps, acc = sqrt(.Machine$double.eps), max.it = 100)
#   matches = match(names(control), names(con))
#   con[matches[!is.na(matches)]] = control[!is.na(matches)]
#   if (any(is.na(matches))) warning("Some names in control not matched: ", paste(names(control[is.na(matches)]), sep = ", "))
# 
#   ## set up output tables with uniform distribution
#   tab0 <- lapply(tab, function(x) array(sum(x)/length(x), dim=dim(x)))
# 
#   ### get neighbours of each clique
#   nbs <- list()
#   for (i in seq_along(clqs)) {
#     tmp <- sapply(seq_along(clqs), function(y) length(intersect(clqs[[i]], clqs[[y]])) > 0)
#     tmp[i] <- FALSE
#     nbs[[i]] = which(tmp)
#   }
#   
#   ## set up main loop
#   it <- 0
#   err <- 1
#   
#   while (it < con$max.it) {
#     it <- it + 1
#     
#     ## go through cliques, matching 
#     for (i in seq_along(clqs)) {
#       tmp <- tab[[i]] * log(tab[[i]]/tab0[[i]])
#       tmp[is.na(tmp)] <- 0
#       KL <- sum(tmp)
#       if (verbose) cat("iteration: ", it, ", clique: ", i, "discrepancy:", signif(KL, 3),"\n")
#       
#       tab0[[i]] <- tab[[i]]
#       for (k in seq_along(nbs[[i]])) {
#         int <- intersect(clqs[[i]], clqs[[nbs[[i]][k]]])
#         wh <- match(int,
#                     clqs[[nbs[[i]][k]]])
#         subtab <- margin.table(tab[[i]], match(int, clqs[[i]]))
#         subtab2 <- margin.table(tab0[[nbs[[i]][k]]], wh)
#         subtab_r <- subtab/subtab2
#         subtab_r[is.na(subtab_r)] <- 0
#         tab0[[nbs[[i]][k]]] <- tab0[[nbs[[i]][k]]]*c(subtab_r[patternRepeat0(which=wh, n = dim(tab0[[nbs[[i]][k]]]))])
#       }
#     }
#     if (verbose) cat("Iteration ", it, " complete\n")
#     
#     ## check for any discrepancy
#     dis <- max(mapply(function(x, y) sqrt(sum(x-y)^2), tab[[i]], tab0[[i]]))
#     if (dis < con$acc) {
#       err <- 0
#       if (verbose) cat("Convergence reached, exiting")
#       break
#     }
#   }
#   
#   ## return relevant sufficient statistics
#   list(tab=tab0, its=it, err=err)
# }
# 
