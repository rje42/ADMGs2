##' Find heads and tails
##' 
##' Lists the (recursive) heads and tails for a CADMG
##' 
##' @param graph object of class \code{graph}, should be a CADMG
##' @param r logical indicating if recursive heads should be used
##' @param byDist logical indicating if results should be returned separated by district
##' @param instrinsic optionally, a list of the relevant intrinsic sets of the ADMG
##' 
##' @details The result is a list containing elements of the 
##' same length named \code{heads}, \code{tails} and (optionally) 
##' \code{intrinsic} giving the heads, tails and intrinsic sets.
##' If \code{byDist} is \code{TRUE} then result is a list 
##' each of which is of the above format, corresponding to the 
##' separate districts of \code{graph}.
##' 
headsTails = function (graph, r = TRUE, byDist = FALSE, intrinsic)
{
  if (byDist) {
    if (missing(intrinsic)) {
      dists <- districts(graph)
      pas <- lapply(dists, function(x) MixedGraphs::adj(graph, x, etype="directed", dir=-1, inclusive=FALSE))
      return(mapply(function(i,j) headsTails(fix(graph,j,i), r=r, byDist=FALSE), dists, pas, SIMPLIFY=FALSE))
    }
    else {
      dists <- lapply(intrinsic, function(x) unique.default(unlist(x)))
      pas <- lapply(dists, function(x) MixedGraphs::adj(graph, x, etype="directed", dir=-1, inclusive=FALSE))
      return(mapply(function(i,j,k) headsTails(graph[i,j], r=r, byDist=FALSE, k), dists, pas, intrisic, SIMPLIFY=FALSE))
    }
  }
  
  ## only coded to here
  if(missing(intrinsic)) intrinsic <- intrinsicSets(graph, r = r, byDist = byDist)
  
  if (byDist) {
    head.list = tail.list = list()
    for (i in seq_along(intrinsic)) {
      if (r) {
        tail.list[[i]] = lapply(intrinsic[[i]], function(v) pa(graph, v))
        head.list[[i]] <- mapply(setdiff, intrinsic[[i]], tail.list[[i]], SIMPLIFY=FALSE)
      }
      else {
        head.list[[i]] = lapply(intrinsic[[i]], function(v) barren(graph, v))
        tail.list[[i]] = list()
        for (j in seq_along(head.list[[i]]))  tail.list[[i]][[j]] = setdiff(union(intrinsic[[i]][[j]], pa(graph, intrinsic[[i]][[j]])), head.list[[i]][[j]])
      }
    }
  }
  else {
    if (r) {
      tail.list = lapply(intrinsic, function(v) pa(graph, v))
      head.list <- mapply(setdiff, intrinsic, tail.list, SIMPLIFY=FALSE)
    }
    else {
      head.list = lapply(intrinsic, function(v) barren(graph, v))
      tail.list = list()
      for (j in seq_along(head.list)) tail.list[[j]] = setdiff(union(intrinsic[[j]], pa(graph, intrinsic[[j]])), head.list[[j]])
    }
  }
  out <- list(heads = head.list, tails = tail.list, intrinsic = intrinsic)
  out
}

##' Finds the set of all intrinsic sets of a CADMG
##' 
##' @param graph object of class \code{mixedgraph}
##' @param r should recursive (i.e. nested) method be used?
##' @param byDist should results be grouped by district?
#
# Algorithm (for recursive case):
#
# step 1: Generate all districts. step 2: Generate all ancestral sets for each
# district. step 3: If an ancestral set is itself a district, it is an intrinsic
# set.  Append it to the output. step 4: If any ancestral sets have multiple
# districts, recurse with that set, and append all intrinsic sets obtained from
# this recursive step to the output.
#
# Note: since we are generating subgraph, we have to reindex nodes in our sets
# from subgraph indexing to the input graph indexing. This is done each time we
# append to the output.  In step 3, we only need to reindex once from the
# district subgraph to the input graph. In step 4 we need to reindex twice, once
# to reindex from the ancestral set subgraph to the district subgraph, and once
# to reindex from the district subgraph to the input graph.
#
# This procedure isn't "efficient," in the sense of generating (some) duplicate
# intrinsic sets, so we remove duplicates in the end.
#
# If the recurse flag is set to false, generate intrinsic sets without
# recursing, which generates sets for the conditional independence factorization
# of mixed graphs.
#
# If the byDist flag is TRUE, index a list of intrinsic sets by the
# outermost district it occurs in.
#
intrinsicSets = function(graph, r = TRUE, byDist = FALSE){
  
  out <- list()
  districts <- districts(graph[random(graph)])
  if (byDist) return(lapply(districts, function(x) intrinsicSets(graph[x],r,FALSE)))
  
  for(i in seq_along(districts)){
    dist <- districts[[i]]
    
    if(r) {
      ## find ancestral subsets of this district
      district_subgraph <- graph[dist]
      an_sets <- anSets(graph[dist])
      
      ## for each ancestral subset find the districts 
      ## and then recurse
      for(an_set in an_sets){
        an_set_subgraph <- district_subgraph[an_set]
        d <- districts(an_set_subgraph)
        
        if(length(d) == 1) {
          out <- c(out, list(an_set))
        }
        if(length(d) > 1) {
          recursed_sets <- Recall(an_set_subgraph, r, FALSE)
          out <- c(out, recursed_sets)
        }
      }
    } # if (r)
    else {
      subs = powerSet(dist)[-1]
      int = logical(length(subs))
      
      ## FOR EACH SUBSET C OF A DISTRICT, CHECK IF IT'S 'INTRINSIC'
      for (j in seq_along(subs)) {
        ## ang IS G_{an C}, SUBGRAPH FORMED BY ANCESTORS OF subs[[j]]
        ang = graph[anc(graph, subs[[j]])]
        ## DISTRICT OF C in ang.
        subdis = dis(ang, subs[[j]])
        
        ## SET C IS INTRINSIC IF IS MAXIMALLY <-> CONNECTED IN G_{an C}
        if (length(subdis) > length(subs[[j]])) int[j] = FALSE # NOT MAXIMAL
        else if(is.na(setmatch(list(subs[[j]]), districts(ang)))) int[j] = FALSE # NOT CONNECTED
        else int[j] = TRUE # OK
      }
      
      out = c(out, subs[int])
    }
  }
  
  # out <- rapply(out, function(x) as.integer(sort.int(x)))
  return(unique(out))
}  # instrinsicSets()



##' Find intrinsic closure of a connected set
intrinsicClosure = function(graph, set, sort=1) {
  #  districts <- districts(graph)
  if (length(set) == 0) return(integer(0))
  if (is(graph, "CADMG") && graph$vtype[set[1]] == "fixed") stop("Not a bidirected connected set of random vertices")
  dist = dis(graph, v=set[1])
  if (!all(set %in% dist)) stop("Not a bidirected connected set of random vertices")
  
  tmp = graph[dist]
  continue = TRUE
  
  while (continue) {
    nv = length(tmp$v)
    tmp = tmp[anc(tmp, set)]
    tmp = tmp[dis(tmp, set)]
    if (length(tmp$v) == nv) break
  }
  out = tmp$v
  
  if (sort > 1) out = sort.int(out)
  return(out)
}


##' Return the partition function for a particular set of vertices
##' 
##' @param graph a CADMG
##' @param ht the output of applying \code{headsTails()} to \code{graph}
##' @param v integer vector of vertices to partition
##' @param r logical indicating whether nested parameterization is being used (only if \code{ht} not provided)
partition = function(graph, ht, v = graph$v, r=TRUE){
  
  if(length(v) == 0) return(list())
  else if (length(v) == 1) return(list(v))
  
  if (missing(ht)) ht = headsTails(graph, r=r)
  headorder = quickSort(ht$intrinsic, subsetOrder)
  
  inc = sapply(ht$heads, function(x) all(x %in% v))
  out = list()
  
  while (any(inc)) {
    wm = max(headorder[inc])
    new = ht$heads[inc & headorder==wm]
    out = c(out, new)
    v = setdiff(v, unlist(new))
    inc = sapply(ht$heads, function(x) all(x %in% v))    
  }
  
  return(out)
}


factorize = function(graph, v = seq_len(n), r = TRUE, ht) {
  n = graph$n
  
  if(length(v) == 0) {
    head.list <- list()
    tail.list <- list()
  }
  else {
    if (missing(ht)) ht = headsTails(graph, r=r)
    head.list <- partition(graph, ht=ht, v, r=r)
    
    wh = setmatch(head.list, ht$heads)
    if (any(is.na(wh))) stop("Error in factorize: some heads not matched")
    tail.list <- ht$tails[wh]
  }
  out = list(heads = head.list, tails = tail.list, vnames=graph$vnames)
  return(out)
}

