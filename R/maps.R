## NOT TESTED AND VERY INEFFICIENT
.getMaps <- 
  function(graph, head_list, tail_list, dist_list, sparse = FALSE, dims = rep(2, n), r = TRUE) {
    
    if (sparse) require(Matrix)
    
    n = length(graph$v)
    # GET HEADS AND TAILS IF NECESSARY
    if (missing(head_list) || missing(tail_list)){
      o <- headsTails(graph, r = r, byDist = TRUE)
      head_list <- lapply(o, function(x) x$heads)
      tail_list <- lapply(o, function(x) x$tails)
    }
    else if(!is.list(head_list[[1]])) stop("Heads must be be by district")
    if (missing(dist_list)) dist_list = districts(graph)
    
    # get heads and tails not by district to use with factorize()
    ht = list(heads=unlist(head_list, recursive=FALSE),
              tails=unlist(tail_list, recursive=FALSE),
              intrinsic=unlist(lapply(o, function(x) x$intrinsic), recursive=FALSE))
    
    # GET UNIONS OF DISTRICTS AND THEIR PARENTS
    pa.dist_list = list()
    for (d in seq_along(dist_list)) pa.dist_list[[d]] = sort(union(dist_list[[d]], pa(graph, dist_list[[d]])))
    
    # PROBMAP A LIST OF MATRICES P_d FOR DISTRICT d.  PARAMAP SAME FOR M_d.
    probmap = paramap = update.list = vector(mode = "list", length = length(dist_list))
    update.list = vector(mode = "list", length = n)
    
    # GET A MAP FOR EACH DISTRICT
    for (d in seq_along(dist_list)) {
      # NUMBER OF TIMES WE NEED q_H|T IS NUMBER OF HEAD COMBINATIONS TIMES TAIL COMBINATIONS INSIDE DISTRICT
      vs = sort(dist_list[[d]])
      patail = sort(setdiff(pa(graph, vs), vs))
      d.subs = combinations(dims[vs])
      all.subs = combinations(dims[sort(c(vs,patail))])
      
      # POSITIONS OF PARAMETERS
      tmp = sapply(head_list[[d]], function(x) prod(dims[x]-1)) * sapply(tail_list[[d]], function(x) prod(dims[x]))
      t = c(1,cumsum(tmp)+1) # STARTING POSITIONS OF q_H|T
      p = tail(t,1)-1 # NUMBER OF PARAMETERS
      
      # MATRICES M, P
      probmap[[d]] = matrix(NA, nrow=prod(dims[c(vs, patail)]), ncol=0)
      paramap[[d]] = matrix(NA, nrow=0, ncol=p)
      
      for (i in seq_len(nrow(d.subs))) {
        fctr = factorize(graph, vs[d.subs[i,] < dims[vs]-1], r=r, ht=ht)
        #      if (length(fctr$heads) == 0) next;
        if (length(fctr$tails) > 0) b.tail = sort.int(unique.default(unlist(fctr$tails)))
        else b.tail = integer(0)
        bt.states = combinations(dims[b.tail])
        tmp1 = matrix(0, nrow=nrow(probmap[[d]]), ncol=prod(dims[b.tail]))
        tmp2 = matrix(0, nrow=ncol(tmp1), ncol=ncol(paramap[[d]]))
        
        # FILL IN NEW M BIT
        # rows IS ROWS WHERE CURRENT TERM IS NEEDED (WITH SIGN)
        rows = apply(all.subs[, match(vs, sort(c(vs,patail))), drop=FALSE], 1, function(x) all((x == d.subs[i,]) | (x == dims[vs]-1)))
        rows = rows*apply(all.subs[, match(vs[d.subs[i,] < dims[vs]-1], sort(c(vs,patail))), drop=FALSE], 1, function(x) (-1)^(sum((x == dims[vs[d.subs[i,] < dims[vs]-1]]-1))))
        
        which.bt = match(b.tail, sort(c(vs,patail)))
        # NOW ADD 1 WHEN TAIL PATTERN MATCHES THAT OF TERM
        for (j in seq_len(ncol(tmp1))) tmp1[apply(all.subs[, which.bt, drop=FALSE], 1, function(x) all(x == bt.states[j,])), j] = 1
        # GET CORRECT SIGN
        tmp1 = tmp1*rows
        
        # P BIT
        for (h in seq_along(fctr$heads)) {   # EACH HEAD IN THE TERM
          which.hd = match(fctr$heads[h], head_list[[d]])
          start = t[which.hd] + sum(c(1, cumprod(dims[fctr$heads[[h]]]-1))*c(d.subs[i, match(fctr$heads[[h]], sort(c(vs)))], 0))
          inc = prod(dims[fctr$heads[[h]]]-1)
          
          # TAIL FOR THIS SPECIFIC HEAD
          this.tail = match(fctr$tails[[h]], b.tail)
          #        cp.tt = cumprod(this.tail)
          cp.tt = cumprod(dims[fctr$tails[[h]]])
          # FOR EACH TAIL STATE, SHOW WHERE THIS GEN. MOEBIUS PARAM. GOES
          for (j in seq_len(nrow(bt.states))) tmp2[j, start+inc*sum(c(1, cp.tt)*c(bt.states[j, this.tail],0))] = 1
        }
        
        probmap[[d]] = cbind(probmap[[d]], tmp1)
        paramap[[d]] = rbind(paramap[[d]], tmp2)
      }
      
      # IF ANY COLUMN OF M IS UNUSED, GET RID OF IT
      tmp = (colSums(abs(probmap[[d]])) == 0)
      probmap[[d]] = probmap[[d]][,!tmp,drop=FALSE]
      paramap[[d]] = paramap[[d]][!tmp,,drop=FALSE]
      
      # LIST OF PARAMETERS TO BE UPDATED AT EACH STAGE
      for (i in dist_list[[d]]) {
        # FOR EACH VERTEX i IN DISTRICT d, WHICH HEADS IS IT IN?
        tmp = which(sapply(head_list[[d]], function(x) i %in% x))
        for (j in tmp) {
          # ALL PARAMETERS ASSOCIATED WITH THESE HEADS MUST BE UPDATED WITH i
          update.list[[i]] = c(update.list[[i]], t[j]:(t[j+1]-1))
        }
      }
      if (sparse) {
        probmap[[d]] = Matrix(probmap[[d]])
        paramap[[d]] = Matrix(paramap[[d]])
      }
    }
    
    list(M = probmap, P = paramap, update.list = update.list, dists = dist_list, pa.dists = pa.dist_list, graph=graph, dim=dims, r=r)
  }


#' Evaluate generalized Moebius parameters
#' 
#' Given a full probability distribution over vertices of a graph, returns the
#' value of the associated generalized Moebius parameters.
#' 
#' 
#' @aliases getMoebius print.mparam
#' @param graph An object of class \code{mixedgraph}, an ADMG.
#' @param ptable An array containing a joint probability distribution over the
#' vertices of \code{graph}.
#' @param r logical - should recursive factorizations be used?
#' @return An object of class \code{mparams}.  That is, a list containing:
#' \item{q}{The Moebius parameters.  \code{q[[d]][[i]][[j]]} is a numeric array
#' containing values of \eqn{P(X_H = x_H \,|\, X_T = j)} for each
#' \eqn{x_H}{x_H}, where \eqn{H}{H} is the \code{i}th head in district
#' \code{d}, and \code{j} indexes the tail states \eqn{x_T}{x_T}.}
#' \item{heads}{List of lists, sorted by district, and each sub-list containing
#' integer vectors of the heads in that district.} \item{tails}{List of lists,
#' containing integer vectors of the tails corresponding to the heads above.}
#' \item{vnames}{Vector of names of vertices in graph.} \item{dim}{Number of
#' categories for each variable.} \item{r}{as input}
#' 
#' @section Warning : This function does not verify that the given distribution
#' satisfies the conditions of the model for the maps being calculated, and
#' thus the distribution will not necessarily be recovered by mapping back
#' using \code{getProbs}.
#' 
#' Note this function will not generally return correct values for parameters
#' in the recursive parametrization unless all probabilities are equal.
#' @author Robin Evans
#' @seealso \code{\link{getProbs}}.
#' @references Evans and Richardson (2010)
#' @examples
#' 
#' data(gr2)
#' 
#' # Distribution of complete independence
#' ptable = array(1/32, rep(2,5))
#' getMoebius(gr2, ptable, r=TRUE)
#' 
getMoebius <-
  function (graph, ptable, r=TRUE) {
    ht = headsTails(graph, r = r, byDist = TRUE)
    head_list <- lapply(ht, function(x) x$heads)
    tail_list <- lapply(ht, function(x) x$tails)
    
    # head_list = tmp$heads; tail_list = tmp$tails
    dims = dim(ptable)
    
    q = vector(mode="list", length=length(ht))
    # EACH DISTRICT i
    for (i in seq_along(q)) {
      q[[i]] = vector(mode="list", length=length(head_list[[i]]))
      # EACH HEAD j
      for (j in seq_along(q[[i]])) {
        # ALL POSSIBLE TAIL ASSIGNMENTS
        vals = combinations(dims[tail_list[[i]][[j]]])
        q[[i]][[j]] = vector(mode="list", length=dim(vals)[1])
        
        # EACH TAIL ASSIGNMENT k
        for (k in seq_len(dim(vals)[1])) {
          # GET CONDITIONAL PROBABILITIES
          tmp = condition.table(ptable, head_list[[i]][[j]], tail_list[[i]][[j]], vals[k,]+1)
          # NEED ALL BUT LAST VALUE IN EACH DIMENSION
          if (is.vector(tmp)) q[[i]][[j]][[k]] = tmp[-length(tmp)]
          else q[[i]][[j]][[k]] = subarray(tmp, lapply(dim(tmp)-1, seq_len))
        }
      }
    }
    out = list(q = q, heads = head_list, tails = tail_list, vnames = graph$vnames, dim=dims, r=r)
    class(out) = "mparam"
    out
  }


#' Return Moebius parameter name strings
#' 
#' Gives a vector of strings corresponding to the meaning of each generalized
#' Moebius parameter: e.g. "p(a, c = 0 | b = 1, d = 0, e = 0)".
#' 
#' The logical \code{blanks} provides for adding whitespace in place of some
#' repeated vertices in the output.
#' 
#' @param mparam An object of class \code{mparam}.
#' @param blanks logical indicating whether to omit repeated names.
#' @return A string valued vector containing the names in head order.
#' @note Note that if the Moebius parametrization is with respect to the
#' recursive factorization, some of the parameters may be reweighted
#' probabilities of the joint distribution.
#' @author Robin Evans
#' @references Evans, R.J. and Richardson, T.S. (2010) - Fitting acyclic
#' directed mixed graphs to binary data. \emph{UAI-10}.
#' @examples
#' 
#' data(gr2)
#' mparams = getMoebius(gr2, makeptable(gr2, values=rep(1/32, 32)))
#' 
#' getMparamsNames(mparams)
#' 
getMparamsNames <-
  function(mparam, blanks = FALSE) {
    vnames = mparam$vnames
    heads = mparam$heads
    tails = mparam$tails
    dims = mparam$dim
    q = mparam$q
    out = character(0)
    
    # FOR EACH HEAD j (IN DISTRICT i)
    for (i in seq_along(q)) for (j in seq_along(q[[i]])) {
      # ALL POSSIBLE TAIL ASSIGNMENTS
      h.len = length(heads[[i]][[j]])
      t.len = length(tails[[i]][[j]])
      
      t.vals = combinations(dims[tails[[i]][[j]]])
      h.vals = combinations(dims[heads[[i]][[j]]]-1)
      
      for (k1 in seq_len(dim(t.vals)[1])) for (k2 in seq_len(dim(h.vals)[1]))  {
        #cat("head =",heads[[i]][[j]], "tail =",tails[[i]][[j]],"\n")
        tmp = "p("
        
        # GET HEAD NAMES AND VALUES (IF-ELSE PRINTS BLANKS IF SAME AS LINE ABOVE AND blanks == TRUE)
        tmp2 = paste(mapply(paste,
                            if(isTRUE(k2 == 1) || !blanks) vnames[heads[[i]][[j]]] else sapply(nchar(vnames[heads[[i]][[j]]]), function(x) paste(rep(" ",x), collapse="")),
                            " = ", h.vals[k2,], c(rep(", ", h.len-1),""), sep=""), collapse="", sep="")
        
        tmp = paste(tmp, tmp2, sep="")
        
        # IF ANYTHING IN TAIL, ADD IN
        if (isTRUE(t.len > 0)) {
          tmp = paste(tmp," | ",sep="")
          tmp2 = paste(mapply(paste, vnames[tails[[i]][[j]]],
                              " = ", t.vals[k1,], c(rep(", ", t.len-1),"") , sep=""), collapse="", sep="")
          tmp = paste(tmp, tmp2, sep="")
        }
        
        tmp = paste(tmp, ")", sep="")
        out = c(out, tmp)
      }
    }
    out
  }

print.mparam <-
  function(x, blanks=FALSE, ...) {
    names = getMparamsNames(x, blanks=blanks)
    values = unlist(x$q)
    
    cat(ifelse(x$r, "Recursive", "Non-recursive"), "ADMG parametrization\n", sep=" ")
    
    for (i in seq_along(values)) {
      cat(names[i], " = ", values[i], "\n", sep="")
    }
  }

#' Constructs probability array for graph
#' 
#' Put joint probabilities in an array with labelled dimensions.
#' 
#' @param graph An ADMG, object of class \code{mixedgraph}.
#' @param dims An integer vector with same length as \code{graph$v} containing the number
#' of categories take by each random variable.
#' @param values A vector of joint probabilities of length \code{prod(dims)}.
#' @return An array containing the probabilities with
#' appropriately labelled dimensions.
#' @author Robin Evans
#' @keywords array
#' @examples
#' 
#' data(gr1)
#' makeptable(gr1)
#' 
makeptable <-
  function(graph, dims = rep(2, n), values = 1/prod(dims)) {
    n <- length(graph$v)
    
    out = array(values, dim=dims)
    if (abs(sum(out)-1) > 1e-10) warning("Probabilities do not sum to 1")
    
    dnames = lapply(dims, function(x) seq_len(x)-1)
    names(dnames) = graph$vnames[graph$v]
    dimnames(out) <- dnames
    
    out
  }

