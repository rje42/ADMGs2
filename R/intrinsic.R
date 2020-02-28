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
##' same length named \code{heads}, \code{tails} and 
##' \code{intrinsic} giving the heads, tails and intrinsic sets.
##' If \code{byDist} is \code{TRUE} then result is a list 
##' each of which is of the above format, corresponding to the 
##' separate districts of \code{graph}.
##' 
##' @export headsTails
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
  
  if (missing(intrinsic)) intrinsic <- intrinsicSets(graph, r = r)
  
  if (r) {
    tail.list = lapply(intrinsic, function(v) pa(graph, v))
    head.list <- mapply(setdiff, intrinsic, tail.list, SIMPLIFY=FALSE)
  }
  else {
    head.list = lapply(intrinsic, function(v) barren(graph, v))
    tail.list = mapply(function(x, y) setdiff(union(x, pa(graph, x)), y), intrinsic, head.list, SIMPLIFY = FALSE)
    # for (j in seq_along(head.list)) tail.list[[j]] = setdiff(union(intrinsic[[j]], pa(graph, intrinsic[[j]])), head.list[[j]])
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
##' List Intrinsic Sets of Summary Graph or ADMG
##' 
##' List Intrinsic Sets of a Summary Graph or Acyclic Directed Mixed Graph, possibly by district.
##' 
##' @param graph object of class \code{mixedgraph}, must be an ADMG.
##' @param r logical, should recursive head definition be used? Defaults to \code{TRUE}
##' @param by_district logical, should intrinsic sets be grouped by district? Defaults to \code{FALSE}
##' @param sort should output be sorted/unique?
##' @param recall logical: is this a recalling of the function (internal use only)
##' 
##' @details \code{intrinsicSets} returns a list of integer vectors, 
##' each being an intrinsic set (or if \code{by_district = TRUE} 
##' a list of lists, each containing the intrinsic sets in a single 
##' district). If \code{r = FALSE} the intrinsic sets are the heads 
##' and their dis-tails.
##' 
##' \code{intrinsicClosure} returns an integer vector containing the closure of set.
##' 
##' @export intrinsicSets
intrinsicSets <- function(graph, r = TRUE, by_district = FALSE, sort=2, recall=FALSE) {
  
  ## on first run, check graph is a summary graph and remove undirected part
  if(!recall) {
    if(!is.SG(graph)) stop("Graph appears not to be a summary graph or ADMG")
    
    un_g <- un(graph)
    if (length(un_g) > 0) {
      clq <- cliques(graph[un_g])
      graph <- graph[-un_g]
    }
  }
  if (recall || length(un_g) == 0) clq <- list()
  
  out <- list()
  districts <- districts(graph[random(graph)])
  if (by_district) return(lapply(districts, function(x) intrinsicSets(graph[x], r, by_district=FALSE, sort=sort, recall=TRUE)))
  
  districts <- districts(graph[random(graph)])
  out <- list()
  
  ## go along each district and find intrinsic sets
  for(i in seq_along(districts)){
    if(by_district) {
      out[[i]] <- list()
    }
    dis <- districts[[i]]
    
    if(r) {
      # if (length(dis) <= 1) {
      if(by_district) {
        out[[i]] <- c(out[[i]], list(dis))
      } 
      else {
        out <- c(out, list(dis))
      }
      
      if (length(dis) <= 1) next
      # }
      
      dis_subgraph <- graph[dis]
      ster <- sterile(dis_subgraph)
      
      for(s in ster) {
        an_set <- setdiff(dis, s)
        
        an_set_subgraph <- dis_subgraph[an_set]
        d <- districts(an_set_subgraph)
        
        recursed_sets <- Recall(an_set_subgraph, r, FALSE, recall=TRUE)
        
        ## add in new sets
        if(by_district){
          out[[i]] <- c(out[[i]], recursed_sets)
          if (length(d) == 1) out[[i]] <- c(out[[i]], list(an_set))
        } else {
          out <- c(out, recursed_sets)
          if (length(d) == 1) out <- c(out, list(an_set))
        }
      }
    } # if (r)
    else {
      subs = powerSet(dis)[-1]
      int = logical(length(subs))
      
      # for each subset C of a district, check if it's intrinsic
      for (j in seq_along(subs)) {
        # ang is G_{an C}, subgraph formed by ancestors of subs[[j]]
        ang = graph[anc(graph, subs[[j]])]
        # district of C in ang.
        subdis = dis(ang, subs[[j]])
        
        # set C is intrinsic if it is <-> connected in G_{an C}
        if (length(subdis) > length(subs[[j]])) int[j] = FALSE # not maximal
        else if(!all(subs[[j]] %in% dis(ang, subs[[j]][1]))) int[j] = FALSE # not connected
        else int[j] = TRUE # OK
      }
      
      if(by_district){
        out[[i]] <- subs[int]
      }
      else out = c(out, subs[int])
    }
  }
  
  ## clean up and finish
  out <- c(clq, out)
    
  if (sort > 0) {
    out <- lapply(out, sort.int)
    out <- unique.default(out)
  }
  out
}

# intrinsicSets = function(graph, r = TRUE, byDist = FALSE) {
#   
#   out <- list()
#   districts <- districts(graph[random(graph)])
#   if (byDist) return(lapply(districts, function(x) intrinsicSets(graph[x],r,FALSE)))
#   
#   for(i in seq_along(districts)){
#     dist <- districts[[i]]
#     
#     if(r) {
#       ## find ancestral subsets of this district
#       district_subgraph <- graph[dist]
#       an_sets <- anSets(graph[dist])
#       
#       ## for each ancestral subset find the districts 
#       ## and then recurse
#       for(an_set in an_sets){
#         an_set_subgraph <- district_subgraph[an_set]
#         d <- districts(an_set_subgraph)
#         
#         if(length(d) == 1) {
#           out <- c(out, list(an_set))
#         }
#         if(length(d) > 1) {
#           recursed_sets <- Recall(an_set_subgraph, r, FALSE)
#           out <- c(out, recursed_sets)
#         }
#       }
#     } # if (r)
#     else {
#       subs = powerSet(dist)[-1]
#       int = logical(length(subs))
#       
#       ## FOR EACH SUBSET C OF A DISTRICT, CHECK IF IT'S 'INTRINSIC'
#       for (j in seq_along(subs)) {
#         ## ang IS G_{an C}, SUBGRAPH FORMED BY ANCESTORS OF subs[[j]]
#         ang = graph[anc(graph, subs[[j]])]
#         ## DISTRICT OF C in ang.
#         subdis = dis(ang, subs[[j]])
#         
#         ## SET C IS INTRINSIC IF IS MAXIMALLY <-> CONNECTED IN G_{an C}
#         if (length(subdis) > length(subs[[j]])) int[j] = FALSE # NOT MAXIMAL
#         else if(is.na(setmatch(list(subs[[j]]), districts(ang)))) int[j] = FALSE # NOT CONNECTED
#         else int[j] = TRUE # OK
#       }
#       
#       out = c(out, subs[int])
#     }
#   }
#   
#   # out <- rapply(out, function(x) as.integer(sort.int(x)))
#   return(unique(out))
# }  # instrinsicSets()


##' Find intrinsic closure of a connected set
##' 
##' @param graph object of class \code{mixedgraph}
##' @param set set of vertices to find intrisic closure ot
##' @param sort 
##' 
##' @export intrinsicClosure
intrinsicClosure = function(graph, set, sort=1) {
  #  districts <- districts(graph)
  if (length(set) == 0) return(integer(0))
  if (is(graph, "CADMG") && any(graph$vtype[set] == "fixed")) stop("Not a bidirected connected set of random vertices")
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


##' @describeIn factorize Return the partition function for a particular set of vertices
## @export partition
partition = function(graph, v = graph$v, r=TRUE, ht, head_order){
  
  if(length(v) == 0) return(list())
  else if (length(v) == 1) return(list(v))
  
  if (missing(ht)) ht = headsTails(graph, r=r)
  if (missing(head_order)) head_order = quickSort(ht$intrinsic, subsetOrder)
  
  inc = sapply(ht$heads, function(x) all(x %in% v))
  out = list()
  
  while (any(inc)) {
    wm = max(head_order[inc])
    new = ht$heads[inc & head_order==wm]
    out = c(out, new)
    v = setdiff(v, unlist(new))
    inc = sapply(ht$heads, function(x) all(x %in% v))    
  }
  
  return(out)
}


##' Factorize a set of nodes into heads and tails
##' 
##' @param graph a CADMG
##' @param v integer vector of vertices to partition
##' @param r logical indicating whether nested parameterization is being used (used only if \code{ht} not provided)
##' @param ht the output of applying \code{headsTails()} to \code{graph}
##' @param head_order optional vector that represents partial order of heads
##' 
## @export factorize
factorize = function(graph, v = graph$v, r = TRUE, ht, head_order) {

  if(length(v) == 0) {
    head.list <- list()
    tail.list <- list()
  }
  else {
    if (missing(ht)) ht = headsTails(graph, r=r)
    if (missing(head_order)) head_order = quickSort(ht$intrinsic, subsetOrder)
    head.list <- partition(graph, ht=ht, v, r=r, head_order=head_order)
    
    wh = setmatch(head.list, ht$heads)
    if (any(is.na(wh))) stop("Error in factorize: some heads not matched")
    tail.list <- ht$tails[wh]
  }
  out = list(heads = head.list, tails = tail.list, vnames=graph$vnames)
  return(out)
}

