## one iteration of fitting algorithm
.doone <- function(dat, graph, params, maps, tol=sqrt(.Machine$double.eps), quietly = TRUE) {
  n = max(graph$v)
  tot = sum(dat[[1]])
  
  if (!quietly) cat("v = ")
  
  for (v in seq_len(n)) {
    if (!quietly) cat(v, " ", sep="")
    d = base::which(sapply(maps$dists, function(x) v %in% x))
    M = maps$M[[d]]; P = maps$P[[d]]; update = maps$update.list[[v]]
    A = .getA(maps, params$q, v, dist=d)
    q = unlist(params$q[[d]])
    
    repeat {
      ll = .logLik2(dat[[d]], q, M, P)
      g = .DlogLik2(dat[[d]], q, M, P, update)
      dir = -g
      dir2 = rep(0,length(q))
      dir2[update] = -g
      
      FUN = function(x) .logLik2(dat[[d]], x, M, P)
      adj.start = 1/max(1,dir2/q)*(1-.Machine$double.eps)
      adj.start = min(adj.start, (1-q[dir2 < -1e-12])/abs(dir2[dir2 < -1e-12]))
      #      adj.start = min(adj.start, (q[dir2 > 1e-12])/dir2[dir2 > 1e-12])
      
      tmp2 = armijo(fun = FUN, x=q, maximise=TRUE, grad=-dir2, adj.start=1/max(1,dir2/q), searchup=FALSE)
      
      diff.ll = tmp2$best - ll
      q = q + tmp2$move
      
      # STOP PROCEDURE IF DIFFERENCE IN LIKELIHOOD TOO SMALL OR MOVES
      # TOO SMALL
      
      finish = (abs(diff.ll/ll) < sqrt(.Machine$double.eps)) ||
        (max(abs(tmp2$move)) < tol)
      if (finish) break
    }
    tmp = 1
    # UPDATE PARAMETER VALUES
    for (i in seq_along(params$q[[d]])) for (j in seq_along(params$q[[d]][[i]])) {
      params$q[[d]][[i]][[j]][] = q[seq.int(from=tmp, length.out=length(params$q[[d]][[i]][[j]]))]
      tmp = tmp + length(params$q[[d]][[i]][[j]])
    }
  }
  
  if (!quietly) cat("\n")
  
  return(params)
}

## Fisher information matrix for parameter fit
.fisher_mixed_fit <- function (fit) {
  q = fit$params$q
  k = length(fit$maps$M)
  # n = fit$graph$n
  p = length(unlist(q))
  dim = fit$dim
  maps = fit$maps

  if (is.array(fit$dat)) freq = c(fit$dat)
  else freq <- fit$dat[,ncol(fit$dat)]
  n_obs = sum(freq)
  probs = rep(1, prod(dim))

  L = matrix(0, nrow = prod(dim), ncol = 0)
  for (i in seq_len(k)) {
    un_q <- unlist(q[[i]])
    probs.i = as.vector(maps$M[[i]] %*% exp(maps$P[[i]] %*%
                                              log(un_q)))
    probs.i = patternRepeat(probs.i, maps$pa.dists[[i]], fit$dim)
    probs = probs * probs.i
    tmp = matrix(as.vector(maps$M[[i]] %*% diag(as.vector(exp(maps$P[[i]] %*%
                                                                log(un_q)))) %*% maps$P[[i]] %*% diag(1/un_q, nrow=length(un_q))),
                 ncol = length(un_q))
    tmp = tmp[patternRepeat(seq_len(prod(fit$dim[maps$pa.dists[[i]]])), maps$pa.dists[[i]], fit$dim), ]
    L = cbind(L, tmp/probs.i)
  }
  L = L * probs # scale
  varn = n_obs * (diag(1/probs, nrow=length(probs)) - 1)
  FIM = t(L) %*% varn %*% L
  FIM
}

.getA <-
  function (maps, q, v, dist) {
    if(missing(dist)) dist = which(sapply(maps$dists, function(x) v %in% x))
    M = maps$M[[dist]]
    P = maps$P[[dist]]
    update = maps$update.list[[v]]
    
    A = matrix(NA, nrow = nrow(M), ncol = length(update) + 1)
    P1 = P[, -update, drop = FALSE]
    P2 = P[, update, drop = FALSE]
    
    # PARAMETERS BEING FIXED
    qq = unlist(q[[dist]])
    qq = qq[-update]
    
    # A[,1] IS 'b' VECTOR FROM PAPER
    if (is(P2, "Matrix")) {
      whR <- (Matrix::rowSums(P2) == 0)
      
      ## For sparse matrices.  Calling slots directly is a bit naughty,
      ## but the alternative is as.vector(), which is very slow.
      A[, 1] = (M[,whR,drop=FALSE] %*% exp(P1[whR,,drop=FALSE] %*% log(qq)))@x
      for (i in seq_along(update)) {
        A[, i+1] = (M[,P[,update[i]]==1,drop=FALSE] %*% exp(P[P[,update[i]]==1,-update,drop=FALSE] %*% log(qq)))@x
      }
    }
    else {
      whR <- (rowSums(P2)==0)
      A[, 1] = M[,whR,drop=FALSE] %*% exp(P1[whR,,drop=FALSE] %*% log(qq))
      for (i in seq_along(update)) {
        A[, i+1] = M[,P[,update[i]]==1,drop=FALSE] %*% exp(P[P[,update[i]]==1,-update,drop=FALSE] %*% log(qq))
      }
    }
    
    A  #* c(const)
  }

# .getMaps <-
#   function(graph, head.list, tail.list, dists.list, sparse = FALSE, dims = rep(2, n), r = TRUE) {
#     
#     if (sparse) requireNamespace(Matrix)
#     
#     n = graph$n
#     # GET HEADS AND TAILS IF NECESSARY
#     if (missing(head.list) || missing(tail.list)){
#       o <- headsTails(graph, r = r, by_district = TRUE)
#       head.list <- o$heads
#       tail.list <- o$tails
#     }
#     else if(!is.list(head.list[[1]])) stop("Heads must be be by district")
#     if (missing(dists.list)) dists.list = districts(graph)
#     
#     # get heads and tails not by district to use with factorize()
#     ht = list(heads=unlist(head.list, recursive=FALSE),
#               tails=unlist(tail.list, recursive=FALSE))
#     
#     # GET UNIONS OF DISTRICTS AND THEIR PARENTS
#     pa.dists.list = list()
#     for (d in seq_along(dists.list)) pa.dists.list[[d]] = sort(union(dists.list[[d]], pa(graph, dists.list[[d]])))
#     
#     # PROBMAP A LIST OF MATRICES P_d FOR DISTRICT d.  PARAMAP SAME FOR M_d.
#     probmap = paramap = update.list = vector(mode = "list", length = length(dists.list))
#     update.list = vector(mode = "list", length = n)
#     
#     # GET A MAP FOR EACH DISTRICT
#     for (d in seq_along(dists.list)) {
#       # NUMBER OF TIMES WE NEED q_{H|T} IS NUMBER OF HEAD COMBINATIONS TIMES TAIL COMBINATIONS INSIDE DISTRICT
#       vs = sort(dists.list[[d]])
#       patail = sort(setdiff(pa(graph, vs), vs))
#       d.subs = combinations(dims[vs])
#       all.subs = combinations(dims[sort(c(vs,patail))])
#       
#       # POSITIONS OF PARAMETERS
#       tmp = sapply(head.list[[d]], function(x) prod(dims[x]-1)) * sapply(tail.list[[d]], function(x) prod(dims[x]))
#       t = c(1,cumsum(tmp)+1) # STARTING POSITIONS OF q_H|T
#       p = tail(t,1)-1 # NUMBER OF PARAMETERS
#       
#       # MATRICES M, P
#       probmap[[d]] = matrix(NA, nrow=prod(dims[c(vs, patail)]), ncol=0)
#       paramap[[d]] = matrix(NA, nrow=0, ncol=p)
#       
#       for (i in seq_len(nrow(d.subs))) {
#         fctr = factorize(graph, vs[d.subs[i,] < dims[vs]-1], r=r, ht=ht)
#         #      if (length(fctr$heads) == 0) next;
#         if (length(fctr$tails) > 0) b.tail = sort.int(unique.default(unlist(fctr$tails)))
#         else b.tail = integer(0)
#         bt.states = combinations(dims[b.tail])
#         tmp1 = matrix(0, nrow=nrow(probmap[[d]]), ncol=prod(dims[b.tail]))
#         tmp2 = matrix(0, nrow=ncol(tmp1), ncol=ncol(paramap[[d]]))
#         
#         # FILL IN NEW M BIT
#         # rows IS ROWS WHERE CURRENT TERM IS NEEDED (WITH SIGN)
#         rows = apply(all.subs[, match(vs, sort(c(vs,patail))), drop=FALSE], 1, function(x) all((x == d.subs[i,]) | (x == dims[vs]-1)))
#         rows = rows*apply(all.subs[, match(vs[d.subs[i,] < dims[vs]-1], sort(c(vs,patail))), drop=FALSE], 1, function(x) (-1)^(sum((x == dims[vs[d.subs[i,] < dims[vs]-1]]-1))))
#         
#         which.bt = match(b.tail, sort(c(vs,patail)))
#         # NOW ADD 1 WHEN TAIL PATTERN MATCHES THAT OF TERM
#         for (j in seq_len(ncol(tmp1))) tmp1[apply(all.subs[, which.bt, drop=FALSE], 1, function(x) all(x == bt.states[j,])), j] = 1
#         # GET CORRECT SIGN
#         tmp1 = tmp1*rows
#         
#         # P BIT
#         for (h in seq_along(fctr$heads)) {   # EACH HEAD IN THE TERM
#           which.hd = match(fctr$heads[h], head.list[[d]])
#           start = t[which.hd] + sum(c(1, cumprod(dims[fctr$heads[[h]]]-1))*c(d.subs[i, match(fctr$heads[[h]], sort(c(vs)))], 0))
#           inc = prod(dims[fctr$heads[[h]]]-1)
#           
#           # TAIL FOR THIS SPECIFIC HEAD
#           this.tail = match(fctr$tails[[h]], b.tail)
#           #        cp.tt = cumprod(this.tail)
#           cp.tt = cumprod(dims[fctr$tails[[h]]])
#           # FOR EACH TAIL STATE, SHOW WHERE THIS GEN. MOEBIUS PARAM. GOES
#           for (j in seq_len(nrow(bt.states))) tmp2[j, start+inc*sum(c(1, cp.tt)*c(bt.states[j, this.tail],0))] = 1
#         }
#         
#         probmap[[d]] = cbind(probmap[[d]], tmp1)
#         paramap[[d]] = rbind(paramap[[d]], tmp2)
#       }
#       
#       # IF ANY COLUMN OF M IS UNUSED, GET RID OF IT
#       tmp = (colSums(abs(probmap[[d]])) == 0)
#       probmap[[d]] = probmap[[d]][,!tmp,drop=FALSE]
#       paramap[[d]] = paramap[[d]][!tmp,,drop=FALSE]
#       
#       # LIST OF PARAMETERS TO BE UPDATED AT EACH STAGE
#       for (i in dists.list[[d]]) {
#         # FOR EACH VERTEX i IN DISTRICT d, WHICH HEADS IS IT IN?
#         tmp = which(sapply(head.list[[d]], function(x) i %in% x))
#         for (j in tmp) {
#           # ALL PARAMETERS ASSOCIATED WITH THESE HEADS MUST BE UPDATED WITH i
#           update.list[[i]] = c(update.list[[i]], t[j]:(t[j+1]-1))
#         }
#       }
#       if (sparse) {
#         probmap[[d]] = Matrix(probmap[[d]])
#         paramap[[d]] = Matrix(paramap[[d]])
#       }
#     }
#     
#     list(M = probmap, P = paramap, update.list = update.list, dists = dists.list, pa.dists = pa.dists.list, graph=graph, dim=dims, r=r)
#   }

.logLik <-
  function(dat, q, maps) {
    k = length(maps$M)
    tmp = 0
    
    for (i in seq_len(k)) {
      log_probs = log(maps$M[[i]] %*% exp(maps$P[[i]] %*% log(unlist(q[[i]]))))
      tmp = tmp + sum(log_probs*c(dat[[i]]))
    }
    
    tmp
  }

.logLik2 <-
  function(dat, q, M, P) {
    if (any(q <= 0)) return(NaN)
    
    tmp = M %*% exp(P %*% log(c(q)))
    
    if (any(tmp <= 0)) return(NaN)
    sum(c(dat) * log(tmp))
  }

.DlogLik <-
  function(dat, q, maps) {
    k = length(maps$M)
    
    out = c()
    
    for (i in seq_len(k)) {
      # DERIVATIVE OF log p(q^k) WITH RESPECT TO q^k.
      un_q <- unlist(q[[i]])
      # tmp = diag(1/as.vector(maps$M[[i]] %*% exp(maps$P[[i]] %*% log(un_q)))) %*%
      #   maps$M[[i]] %*% diag(as.vector(exp(maps$P[[i]] %*% log(un_q)))) %*% maps$P[[i]] %*% diag(1/un_q, nrow=length(un_q))
      ePlq <- as.vector(exp(maps$P[[i]] %*% log(un_q)))
      tmp = ((1/c(as.vector(maps$M[[i]] %*% ePlq)))*maps$M[[i]]) %*% 
        (c(ePlq)*maps$P[[i]]) %*% diag(1/un_q, nrow=length(un_q))
      # GET MATRIX TO RIGHT DIMENSION
      tmp = c(dat[[i]]) %*% tmp
      
      out <- c(out, as.vector(tmp))
      # out = as.vector(c(out, tmp[-1]))
    }
    out
  }

.DlogLik2 <-
  function(dat, q, M, P, update) {
    tmp = exp(P %*% log(q))
    tmp2 = (1/as.vector(M %*% tmp) * (M %*% (as.vector(tmp) * (P[,update] * 1/rep(q[update],each=nrow(P))))))
    # tmp2 = (1/as.vector(M %*% tmp) * (M %*% (as.vector(tmp) * (P %*% diag(1/q, nrow=length(1/q))[,update]))))
    out = c(dat) %*% tmp2
    
    as.vector(out)
  }
