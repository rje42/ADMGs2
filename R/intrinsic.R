##' Find heads and tails
##' 
##' Lists the (recursive) heads and tails for a CADMG
##' 
##' @param graph object of class \code{graph}, should be a CADMG
##' @param r logical indicating if recursive heads should be used
##' @param by_district logical indicating if results should be returned separated by district
##' @param sort should output be unique/sorted?
##' @param intrinsic optionally, a list of the relevant intrinsic sets of the ADMG
##' @param max_head optional maximum size for each head
##' 
##' @details The result is a list containing elements of the 
##' same length named \code{heads}, \code{tails} and 
##' \code{intrinsic} giving the heads, tails and intrinsic sets.
##' If \code{by_district} is \code{TRUE} then result is a list 
##' each of which is of the above format, corresponding to the 
##' separate districts of \code{graph}.
##' 
##' @export headsTails
headsTails = function (graph, r = TRUE, by_district = FALSE, sort=1, intrinsic, max_head)
{
  ung <- un(graph)
  if (length(ung) > 0) gr2 <- graph[-ung]
  else gr2 <- graph
  
  ## if no vertices left, return empty lists
  if (nv(gr2) == 0) return(list(heads=list(), tails=list(), intrinsic=list()))
  
  ## otherwise continue
  if (by_district && r) {
      if (missing(intrinsic)) {
        dists <- districts(gr2)
        pas <- lapply(dists, function(x) MixedGraphs::adj(graph, x, etype="directed", dir=-1, inclusive=FALSE))
        out <- mapply(function(i,j) headsTails(fix(graph,j,i), r=r, by_district=FALSE), dists, pas, SIMPLIFY=FALSE)
        
        
        if (missing(max_head)) return(out)
        else {
          for (i in seq_along(dists)) {
            kp <- lengths(out[[i]]$heads) <= max_head
            out[[i]]$heads <- out[[i]]$heads[kp]
            out[[i]]$tails <- out[[i]]$tails[kp]
          }
          return(out)
        }
      }
      else {
        dists <- lapply(intrinsic, function(x) unique.default(unlist(x)))
        pas <- lapply(dists, function(x) MixedGraphs::adj(graph, x, etype="directed", dir=-1, inclusive=FALSE))
        out <- mapply(function(i,j,k) headsTails(graph[i,j], r=r, by_district=FALSE, intrinsic = k), 
                      dists, pas, intrinsic, SIMPLIFY=FALSE)
        
        if (missing(max_head)) return(out)
        else {
          kp <- lengths(out$heads) <= max_head
          out$heads <- out$heads[kp]
          out$tails <- out$tails[kp]
          return(out)
        }
      }
  }
  else if (by_district && !r) {
    if (missing(max_head)) intrinsic <- intrinsicSets(gr2, r=FALSE, by_district = TRUE)
    else intrinsic <- intrinsicSets2(gr2, r=FALSE, by_district = TRUE, maxbarren = max_head)
    
    head.list = tail.list = list()
    for (i in seq_along(intrinsic)) {
      head.list[[i]] = lapply(intrinsic[[i]], function(v) barren(graph, v))
      tail.list[[i]] = list()
      for (j in seq_along(head.list[[i]])) tail.list[[i]][[j]] = setdiff(union(intrinsic[[i]][[j]], pa(graph, intrinsic[[i]][[j]])), head.list[[i]][[j]])
    }
    
    out <- list(heads = head.list, tails = tail.list, intrinsic = intrinsic)
    return(purrr::transpose(out))
  }
#  }
  
  if (missing(intrinsic) && !r) {
    if (missing(max_head)) intrinsic <- intrinsicSets(gr2, r = FALSE, sort=sort)
    else intrinsic <- intrinsicSets2(gr2, r = FALSE, sort=sort, maxbarren = max_head)
  }
  else if (missing(intrinsic)) intrinsic <- intrinsicSets(gr2, r = TRUE, sort=sort)
  
  if (r) {
    tail.list = lapply(intrinsic, function(v) pa(graph, v))
    head.list <- mapply(setdiff, intrinsic, tail.list, SIMPLIFY=FALSE)
  }
  else {
    head.list = lapply(intrinsic, function(v) barren(graph, v))
    tail.list = mapply(function(x, y) setdiff(union(x, pa(graph, x)), y), intrinsic, head.list, SIMPLIFY = FALSE)
    # for (j in seq_along(head.list)) tail.list[[j]] = setdiff(union(intrinsic[[j]], pa(graph, intrinsic[[j]])), head.list[[j]])
  }
  
  if (!missing(max_head)) {
    wh <- lengths(head.list) <= max_head
    head.list <- head.list[wh]
    tail.list <- tail.list[wh]
  }
  
  if (sort > 1) {
    head.list <- lapply(head.list, sort.int)
    tail.list <- lapply(tail.list, sort.int)
    intrinsic <- lapply(intrinsic, sort.int)
  }
  if (sort > 2) {
    ord <- order(sapply(head.list, function(x) sum(2^(x-1))))
    head.list <- head.list[ord]
    tail.list <- tail.list[ord]
    intrinsic <- intrinsic[ord]
  }
  
  out <- list(heads = head.list, tails = tail.list, intrinsic = intrinsic)
  out
}


# 
# ##' Get list of heads and tails
# ##' 
# ##' @param graph a \code{mixedgraph} object
# ##' @param r logical, should recursive heads be used?
# ##' @param by_district logical, should these be grouped by district?
# ##' @param set_list list of intrinsic sets
# ##' @param max_head largest head size to consider (\code{r=FALSE} only)
# ##' 
# ##' @details Returns a list of heads and their corresponding tails.
# ##' 
# ##' @export headsTails
# headsTails <- function (graph, r = TRUE, by_district = FALSE, set_list, max_head) 
# {
#   if(missing(set_list)) {
#     if (missing(max_head)) set_list <- intrinsicSets(graph, r = r, by_district = by_district)
#     else {
#       if (r) stop("Function does not support maximizing the size of recursive heads")
#       set_list <- intrinsicSets2(graph, r = r, by_district = by_district, maxbarren = max_head)
#     }
#   }
#   
#   if (by_district) {
#     head.list = tail.list = list()
#     for (i in seq_along(set_list)) {
#       if (r) {
#         tail.list[[i]] = lapply(set_list[[i]], function(v) pa(graph, v))
#         head.list[[i]] <- mapply(setdiff, set_list[[i]], tail.list[[i]], SIMPLIFY=FALSE)
#       }
#       else {
#         head.list[[i]] = lapply(set_list[[i]], function(v) barren(graph, v))
#         tail.list[[i]] = list()
#         for (j in seq_along(head.list[[i]]))  tail.list[[i]][[j]] = setdiff(union(set_list[[i]][[j]], pa(graph, set_list[[i]][[j]])), head.list[[i]][[j]])
#       }
#     }
#   }
#   else {
#     if (r) {
#       tail.list = lapply(set_list, function(v) pa(graph, v))
#       head.list <- mapply(setdiff, set_list, tail.list, SIMPLIFY=FALSE)
#     }
#     else {
#       head.list = lapply(set_list, function(v) barren(graph, v))
#       tail.list = list()
#       for (j in seq_along(head.list))  tail.list[[j]] = setdiff(union(set_list[[j]], pa(graph, set_list[[j]])), head.list[[j]])
#     }
#   }
#   out <- list(heads = head.list, tails = tail.list)
#   out
# }
# 


## Finds the set of all intrinsic sets of a CADMG
## 
## @param graph object of class \code{mixedgraph}
## @param r should recursive (i.e. nested) method be used?
## @param byDist should results be grouped by district?
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
  if (by_district && r) return(lapply(districts, 
                                 function(x) intrinsicSets(graph[x], r, by_district=FALSE, sort=sort, recall=TRUE)))
  
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
        
        recursed_sets <- Recall(an_set_subgraph, r=r, FALSE, recall=TRUE)
        
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
    } # else [!r]
  }
  
  ## clean up and finish
  out <- c(clq, out)
    
  if (by_district) {
    if (sort > 0) {
      rapply(out, sort.int, how = "replace")
      out <- lapply(out, unique.default)
    }
  }
  else {
    if (sort > 0) {
      out <- lapply(out, sort.int)
      out <- unique.default(out)
    }
    if (sort > 2) {
      ord <- order(sapply(out, function(x) sum(2^(x-1))))
      out <- out[ord]
    }
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


#' @param maxbarren Maximum number of barren nodes (i.e. head size) to consider
#' @describeIn intrinsicSets Alternative method for non-recursive heads only
#' @export intrinsicSets2
intrinsicSets2 <- function(graph, r = TRUE, by_district = FALSE, maxbarren, sort=1) {
  districts <- districts(graph)
  
  if (missing(maxbarren) || maxbarren > nv(graph)) maxbarren = nv(graph)
  out <- list()
  
  if(!is.ADMG(graph)) stop("Graph appears not to be an ADMG") # could extend to MEGs
  
  subs <- anSets2(graph, maxbarren = maxbarren, same_dist = TRUE)
  if (by_district) d <- sapply(subs, function(x) subsetmatch(x[1], districts))
  
  if(r) {
    # stop("function does not work for recursive heads")
    
    out <- lapply(graph$v, function(set) intrinsicClosure(graph, set, r=TRUE))
    n <- length(graph$vnames)
    vnam <- sapply(out, paste, collapse="")
    liv <- rep(TRUE, n)[graph$v]
    liv[lengths(liv) >= maxbarren] <- FALSE
    
    bid <- matrix(0, n, n)
    for (i in seq_len(n)[-1]) for (j in seq_len(i-1)) bid[i,j] <- bid[j,i] <- 1*any(sib(graph, out[[i]]) %in% out[[j]])
    sub <- matrix(0, n, n)
    for (i in seq_len(n)[-1]) for (j in seq_len(i-1)) sub[j,i] <- 1*is.subset(out[[j]],out[[i]])
    
    edges <- list(directed=sub, bidirected=bid)
    
    graph2 <- mixedgraph(v=seq_len(n), edges=edges, vnames=vnam)
    new <- any(bid > 0)
    
    while (new) {
      new <- FALSE
      bidi <- which(graph2$edges$bidirected > 0 & upper.tri(graph2$edges$bidirected))
      for (i in seq_len(nrow(bidi))) {
        if (liv[bidi[i,1]] && liv[bidi[i,2]]) {
          out[[n+1]] <- intrinsicClosure(graph, c(out[[bidi[i,1]]], out[[bidi[i,2]]]))
          graph2 <- MixedGraphs::addNodes(graph2, 1, vnames=paste(out[[n+1]], collapse=""))
          n <- n+1
          
          for (j in seq_len(n-1)) {
            graph2$edges$directed[j,n] <- 1*is.subset(out[[j]],out[[n]])
            graph2$edges$directed[n,j] <- 1*is.subset(out[[n]],out[[j]])
          }
          
          new <- TRUE
        }
        
        if (any(liv[dec(graph2, out[[bidi[i,1]]])]) && 
            any(liv[dec(graph2, out[[bidi[i,2]]])])) {
          
          for (j in dec(graph2, out[[bidi[i,1]]])) for (k in dec(graph2, out[[bidi[i,2]]])) {
            
            
            if (j == 1 && k == 1) next
            
            out[[n+1]] <- intrinsicClosure(graph, c(out[[j]], out[[k]]))
            graph2 <- MixedGraphs::addNodes(graph2, 1, vnames=paste(out[[n+1]], collapse=""))
            n <- n+1
            liv[n] <- TRUE
            
            for (j in seq_len(n-1)) {
              graph2$edges$directed[j,n] <- 1*is.subset(out[[j]], out[[n]])
              graph2$edges$directed[n,j] <- 1*is.subset(out[[n]], out[[j]])
            }
            
            new <- TRUE
            
          }
        }
        
        graph2 <- removeEdges(graph2, list(bidirected=list(bidi[i,])))
      }
    }
    
    stop("function does not work for recursive heads")
    
    ### insert work here
    
  } # if (r)
  else {
    # subs = powerSet(dis, maxbarren)[-1]
    int = logical(length(subs))
    
    # FOR EACH SUBSET C OF A DISTRICT, CHECK IF IT'S 'INTRINSIC'
    for (j in seq_along(subs)) {
      # ang IS G_{an C}, SUBGRAPH FORMED BY ANCESTORS OF subs[[j]]
      bar <- barren(graph, subs[[j]])
      dist <- dis(graph[subs[[j]]], bar[1])
      dist2 <- dis(graph[subs[[j]]], bar)
      if (setequal(dist2, dist)) {
        out <- c(out, list(dist))
        int[j] = TRUE
      }
      else int[j] = FALSE
      
      # 
      # ang = graph[anc(graph, subs[[j]])]
      # # DISTRICT OF C in ang.
      # subdis = dis(ang, subs[[j]])
      # 
      # # SET C IS INTRINSIC IF IS MAXIMALLY <-> CONNECTED IN G_{an C}
      # if (length(subdis) > length(subs[[j]])) int[j] = FALSE # NOT MAXIMAL
      # else if(!all(subs[[j]] %in% dis(ang, subs[[j]][1]))) int[j] = FALSE # NOT CONNECTED
      # else int[j] = TRUE # OK
    }
  }
  
  # out <- subs[int]
  
  if(by_district){
    out <- tapply(out, INDEX = d[int], FUN=list)
    names(out) <- NULL
    
    if (sort > 1) out <- rapply(out, sort.int, how="list")
  } 
  else if (sort > 1) {
    out <- lapply(out, sort.int)
  }
  
  out
}


##' Find intrinsic closure of a connected set
##' 
##' @param graph object of class \code{mixedgraph}
##' @param set set of vertices to find intrisic closure of
##' @param r logical: should recursive heads be used?
##' @param sort if sort > 1 then output is sorted
##' 
##' @export intrinsicClosure
intrinsicClosure = function(graph, set, r=TRUE, sort=1) {
  #  districts <- districts(graph)
  if (length(set) == 0) return(integer(0))
  if (is(graph, "CADMG") && any(graph$vtype[set] == "fixed")) stop("Not a bidirected connected set of random vertices")
  dist = dis(graph, v=set[1])
  if (!all(set %in% dist)) stop("Not a bidirected connected set of random vertices")
  
  ## non-recursive case
  if (!r) {
    ancs <- anc(graph, set)
    S <- dis(graph[ancs], set, sort=sort)
    if (length(districts(graph[S])) <= 1) return(S)
    else stop("Not contained in a single 'intrinsic' set")
  }

  ## deal with recursive case
  tmp2 <- graph$v
  tmp <- graph[dist]
  continue <- TRUE
  
  if (length(tmp$v) > length(set)) {
    while (continue) {
      tmp2 <- tmp$v
      nv = length(tmp2)
      tmp = tmp[dis(tmp, set)]
      tmp = tmp[anc(tmp, set)]
      if (length(tmp$v) == nv) break
    }
  }
  out = tmp$v
  
  if (sort > 1) out = sort.int(out)

  return(out)
}

##' Get C^*
##' 
##' @param graph ADMG object of class \code{mixedgraph}
##' @param C \code{intrinsic set}
##' @param int optional list of intrinsic sets
##' 
##' @details Obtains largest intrinsic set of which the 
##' intrinsic set \code{C} is an ancestral district.  If \code{C} is
##' an ancestral district, then the code just returns 
##' the vertices of \code{graph}.
##' 
##' @export
ancDisClos <- function (graph, C, int) {
  
  ## deal with top level sets
  tmp <- graph[anc(graph, C)]
  tmp <- dis(tmp, C[1])
  if (setequal(C, tmp)) return(graph$v)

  ## otherwise, get a list of the intrinsic sets
  if (missing(int)) int <- intrinsicSets(graph, sort = 3)

  wh <- sapply(int, function(x) is.subset(C, x))
  int <- int[wh]

  int <- int[lengths(int) > length(C)]
  out <- rep(FALSE, length(int))
  
  for (i in seq_along(int)) {
    tmp <- graph[int[[i]]]
    tmp <- tmp[anc(tmp, C)]
    if (setequal(C, dis(tmp, C[1]))) out[i] <- TRUE
    # if (setmatch(C, districts(graph), nomatch = 0) > 0) out[i] <- TRUE
  }

  int <- int[out]
  kp <- list()
  int <- int[order(sapply(int, function(x) sum(2^x)))]
  
  if (length(int) == 1) kp <- int
  else if (length(int) > 1) {
    while (length(int) > 1) {
      kp <- c(kp, last(int))
      int <- int[c(sapply(int[-length(int)], function(x) !is.subset(x, int[[length(int)]])), FALSE)]
    }  
  }  
  else stop("No intrinisic sets containing C")
  
  Cstar <- Reduce(intersect, kp)
  Cstar <- dis(graph[Cstar], C[1])
  
  # Cstar <- Reduce(union, int[out])
  # if (setmatch(list(Cstar), c(int, list(graph$v)), nomatch = 0) == 0) stop("No unique maximal intrinsic set")
  
  # tmp <- graph[anc(graph, C)]
  # tmp <- dis(tmp, C[1])
  # if (setequal(C, tmp)) return(graph$v)
  # if (is.subset(tmp, C)) stop("Not an intrinsic set")
  # 
  # Cstar <- integer(0)
  # if (length(graph$v) == 0) return(integer(0))
  # change <- TRUE
  # 
  # while (change) {
  #   if (setequal(graph$v, Cstar)) stop("Not an intrinsic set")
  #   graph <- graph[dis(graph, C)]
  #   Cstar <- graph$v
  #   graph <- graph[anc(graph, C)]
  #   
  #   if(setequal(C, graph$v)) break
  # }

  Cstar
}

# ##' @describeIn intrinsicSets Get the intrinsic closure of a set
# ##' @param set set to find intrinsic closure of
# ##' @export intrinsicClosure
# intrinsicClosure <- function(graph, set, r=TRUE, sort=1) {
#   
#   if (length(set) == 0) return(integer(0))
#   
#   subgraph <- graph
#   new = TRUE
#   ans <- graph$v
#   
#   if (!r) {
#     ancs <- anc(graph, set)
#     S <- dis(graph[ancs], set, sort=sort)
#     if (length(districts(graph[S])) <= 1) return(S)
#     else stop("Not contained in a single intrinsic set")
#   }
# 
#   while (new) {
#     new = FALSE
#     dists <- districts(subgraph)
#     
#     wh <- subsetmatch(list(set[1]), dists)
#     if (all(set %in% dists[[wh]])) {
#       subgraph <- subgraph[dists[[wh]]]
#     }
#     else stop("Not contained in a single intrinsic set")
#     
#     ans <- anc(subgraph, set)
#       
#     if (setequal(ans, subgraph$v)) break
#     else {
#       subgraph <- subgraph[ans]
#       new = TRUE
#     }
#   }
#   
#   return(subgraph$v)
#   # stop("Error in intrinsicClosure - perhaps 'set' not contained in a single district?")
# }


##' @describeIn factorize Return the partition function for a particular set of vertices
## @export partition
partition = function(graph, v = graph$v, r=TRUE, ht, head_order) {
  
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

##' Gives partition/factorization
##' 
##' Uses order for speed
##' returns integer value of heads from provided list
##' 
##' @param graph object of class \code{mixedgraph}
##' @param heads list of heads 
##' @param v set of vertices to partition
##' @param r logical: should recursive parameterization be used?
##' @param head.order numeric vector in same order as heads
##' 
##' @export partition0
partition0 = function(graph, heads, v = seq_len(graph$n), r=TRUE, head.order) {
  
  wh = rep.int(TRUE, length(heads))
  out = numeric(0)
  
  while (length(v) > 0) {
    wh2 = fsapply(heads[wh], function(x) is.subset(x,v))
    wh[wh] = wh[wh] & wh2
    maxval = max(head.order[wh])
    new = which(wh & head.order==maxval)
    out = c(out, new)
    v = setdiff(v, unlist(heads[new]))
  }
  
  return(out)
}

##' @describeIn partition0 Give factorization of heads and tails
##' @param ht list of heads and tails
##' @export factorize0
factorize0 = function (graph, v = seq_len(n), r = TRUE, ht, head.order) {
  n = graph$n
  if (length(v) == 0) {
    return(numeric(0))
  }
  else {
    out <- partition0(graph, ht$heads, v, r = r, head.order=head.order)
  }
  
  return(out)
}

