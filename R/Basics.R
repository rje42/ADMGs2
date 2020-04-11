graphOptionsEnv <- MixedGraphs:::graphOptionsEnv

## List of vertex types, which can be expanded [not currently supported]
## type     : character name of vertex type
assign("vertexTypesDF", data.frame(type=c("random", "fixed"),
                                   hidden=c(FALSE, FALSE)), 
       envir=graphOptionsEnv)

##' Print method for a CADMG object
##' 
##' @param x a \code{mixedgraph} that is also a \code{CADMG}
##' @param ... other arguments to print
##' 
##' @details Prints a CADMG nicely onto the standard output
##' 
##' @export
print.CADMG = function(x, ...) {
  w = x$v[x$vtypes=="fixed"]
  v = x$v[x$vtypes=="random"]
  nv = length(v)
  nw = length(w)
  cat("CADMG with ", nv, ifelse(nv == 1, " random vertex", " random vertices"),
      ifelse(nv == 0,"",":  "), sep="")
  cat(x$vnames[v], "\n", sep="  ")
  cat("       and ", nw, ifelse(nw == 1, " fixed vertex ", " fixed vertices "),
      ifelse(nw == 0,"",":  "), sep="")
  cat(x$vnames[w], "\n\n", sep="  ")
  
  # get list of edge symbols
  whEdge <- match(names(x$edges), edgeTypes()$type)
  edgeSymb <- edgeTypes()$char
  
  for (i in seq_along(x$edges)) {
    if (is.matrix(x$edges[[i]])) {
      tmp <- cbind(row(x$edges)[x$edges > 0], col(x$edges)[x$edges > 0])
      if(!edgeTypes()$directed[whEdge[i]]) tmp = tmp[tmp[,1] < tmp[,2]]
      
      for (j in seq_len(nrow(tmp))) {
        cat(x$vnames[tmp[j,1]], edgeSymb[whEdge[i]],
            x$vnames[tmp[j,2]], "\n", sep=" ")
      }
      cat("\n")    
    }
    else if (is.data.frame(x$edges[[i]])) {
      for (j in seq_len(nrow(x$edges[[i]]))) {
        cat(x$vnames[x$edges[[i]][j,1]], edgeSymb[whEdge[i]],
                     x$vnames[x$edges[[i]][j,2]], "\n", sep=" ")
      }
      cat("\n")
    }      
    else if (is.list(x$edges[[i]])) {
      if (!is.null(x$edges[[i]]) && length(x$edges[[i]]) > 0) {
        for (j in seq_along(x$edges[[i]])) {
          cat(x$vnames[x$edges[[i]][[j]][1]], edgeSymb[whEdge[i]],
              x$vnames[x$edges[[i]][[j]][2]], "\n", sep=" ")
        }
        cat("\n")
      }
    }
  } 
  
  invisible(x)
}

##' Fix a node in a CADMG
##' 
##' @param graph (C)ADMG object of class mixedgraph
##' @param w set of nodes to fix
##' @param v set of random nodes to retain (defaults to all not fixed)
##'
##' Given a CADMG or ADMG returns the CADMG obtained by 
##' fixing the vertices in \code{w}.
##' 
##' @export fix
fix = function(graph, w, v) {
  if (missing(v)) v = setdiff(graph$v, w)
  else if (length(intersect(w,v))) stop("Fixed and random vertex sets must be disjoint")
  else graph = graph[c(w,v)]
  
  if (any(!(c(w,v) %in% graph$v))) stop(paste("No node ", paste(setdiff(c(w,v), graph$v), sep=", "), " found", sep=""))
  
  graph <- withEdgeList(graph)
  
  e = graph$edges
  if (length(e$directed)) {
    rmv = (matrix(unlist(e$directed), nrow=2)[2,] %in% w)
    e$directed = e$directed[!rmv]
  }
  if (length(e$bidirected)) {
    bied = matrix(unlist(e$bidirected), nrow=2)
    rmvbi = (bied[1,] %in% w) | (bied[2,] %in% w)
    e$bidirected = e$bidirected[!rmvbi]
  }
  if (length(e$undirected)) {
    stop("Only valid for CADMGs")
  }
  vtype = rep("fixed", length(c(v,w)))
  vtype[graph$v %in% v] = "random"
  out <- list(v=graph$v, edges=e, vnames=graph$vnames, vtypes=vtype)
  class(out) <- c("CADMG", "mixedgraph")
  out
}

##' Get vector of fixed or random nodes in CADMG
##' 
##' @param graph a \code{mixedgraph} that is also a \code{CADMG}
##' 
##' @export
random = function(graph) {
  if (is(graph, "CADMG")) {
    return(graph$v[(graph$vtypes == "random")])
  }
  else return(graph$v)
}

##' @describeIn random get fixed nodes
##' @export
fixed = function(graph) {
  if (is(graph, "CADMG")) {
    return(graph$v[(graph$vtypes == "fixed")])
  }
  else return(integer(0))
}

# ##' Compute non-recursive dis-tail
# ##' 
# ##' @param graph a \code{mixedgraph} that is also a \code{(C)ADMG}
# ##' @param head set to calculate \code{distail} for
# ##' 
# distail = function(graph, head) {
#   subgraph = graph[anc(graph, head)]
#   vs = dis(subgraph, head)
#   out = setdiff(vs, head)
#   
#   out
# }



# print.mixed.factorization = function (x, ...) {
#   vnames = x$vnames; heads = x$heads; tails = x$tails
#   for (i in seq_along(heads)) {
#     cat("p(")
#     cat(vnames[heads[[i]]],sep=", ")
#     if (length(tails[[i]]) > 0) {
#       cat(" | ")
#       cat(vnames[tails[[i]]],sep=", ")
#     }
#     cat(") ")
#   }
#   if (length(heads) == 0L) print(1)
#   cat("\n")
# }




# 
# gr1 = graphCr("1->3<->2->4<->1")
# gr2 = graphCr("5->3,4->2->1","1<->3<->2<->5<->4<->3")
# grVerma = graphCr("1->2->3->4,2<->4,1->3")
# grIV = graphCr("1->2->3<->2")
# # 
# # gr0 = makeGraph(0)
# 
# gss = structure(list(Trust = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
#                                0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
#                                1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
#                                0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
#                                1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
#                                0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
#                                1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
#                      Helpful = c(0, 0, 1, 1, 0,
#                                  0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
#                                  1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1,
#                                  1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
#                                  0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0,
#                                  0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
#                                  1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
#                      MemUnion = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1,
#                                   1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
#                                   0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
#                                   0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1,
#                                   1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0,
#                                   0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
#                                   1),
#                      MemChurch = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
#                                    1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
#                                    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
#                                    1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
#                                    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
#                                    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
#                                    1, 1, 1, 1, 1, 1, 1, 1),
#                      ConLegis = c(0, 0, 0, 0, 0, 0, 0, 0,
#                                   0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
#                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
#                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                   1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
#                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
#                      ConClerg = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                   1),
#                      ConBus = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                                 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                 1, 1, 1, 1, 1, 1, 1),
#                      freq = c(18, 4, 5, 5, 79, 47, 17, 30, 8,  9, 1, 15, 88,
#                               55, 22, 79, 22, 11, 10, 13, 194, 95, 33, 77, 31,
#                               10, 13, 23, 179, 82, 58, 122, 7, 5, 1, 3, 40, 27,
#                               11, 23, 9, 10, 1, 12, 68, 56, 33, 73, 15, 13, 6,
#                               14, 188, 117, 52, 100, 32, 28, 22, 35, 366, 185,
#                               120, 312, 7, 5, 2, 6, 62, 32, 11, 48, 5, 9, 2, 12,
#                               38, 37, 11, 64, 40, 26, 17, 34, 270, 187, 73, 281,
#                               25, 33, 11, 50, 202, 216, 84, 356, 5, 2, 3, 11,
#                               51, 32, 17, 59, 15, 18, 7, 33, 104, 79, 40, 172,
#                               74, 62, 27, 108, 603, 469, 177, 654, 199, 181, 84,
#                               305, 1002, 920, 460, 1818)),
#                 .Names = c("Trust",
#                            "Helpful", "MemUnion", "MemChurch", "ConLegis",
#                            "ConClerg", "ConBus",
#                            "freq"), row.names = c(NA, -128L), class = "data.frame")
# 
# gss_small = structure(list(Trust = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
#                            Helpful = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
#                            MemUnion = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
#                            MemChurch = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
#                            freq = c(188, 128, 71, 194, 1487, 1006, 391, 1272,
#                                     324, 298, 141, 485, 2047, 1630, 828, 2996)),
#                       .Names = c("Trust", "Helpful", "MemUnion", "MemChurch",
#                                  "freq"), row.names = c(NA, -16L), class = "data.frame")
# 
# twins = expand.grid(list(c(0,1))[rep.int(1,4)])
# twins = cbind(twins, c(288, 80, 92, 51, 15, 9, 7, 10, 8, 4, 8, 9, 3, 2, 4, 7))
# colnames(twins) = c("A1","A2","D1","D2","freq")
# 
