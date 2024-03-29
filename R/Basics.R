graphOptionsEnv <- MixedGraphs:::graphOptionsEnv

## List of vertex types, which can be expanded [not currently supported]
## type     : character name of vertex type
## hidden   : logical for whether vertex is unobserved
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
  w = x$v[na.omit(x$vtypes=="fixed")]
  v = x$v[na.omit(x$vtypes=="random")]
  nv = length(v)
  nw = length(w)
  cat("CADMG with ", nv, ifelse(nv == 1, " random vertex", " random vertices"),
      ifelse(nv == 0,"",":  "), sep="")
  cat(x$vnames[v], "\n", sep="  ")
  cat("       and ", nw, ifelse(nw == 1, " fixed vertex ", " fixed vertices "),
      ifelse(nw == 0,"",":  "), sep="")
  cat(x$vnames[w], "\n\n", sep="  ")

  MixedGraphs:::print.edgeList(x$edges, x$vnames)
  #   
  # # get list of edge symbols
  # whEdge <- match(names(x$edges), edgeTypes()$type)
  # edgeSymb <- edgeTypes()$char
  # 
  # for (i in seq_along(x$edges)) {
  #   if (is.matrix(x$edges[[i]])) {
  #     tmp <- cbind(row(x$edges)[x$edges > 0], col(x$edges)[x$edges > 0])
  #     if(!edgeTypes()$directed[whEdge[i]]) tmp = tmp[tmp[,1] < tmp[,2]]
  #     
  #     for (j in seq_len(nrow(tmp))) {
  #       cat(x$vnames[tmp[j,1]], edgeSymb[whEdge[i]],
  #           x$vnames[tmp[j,2]], "\n", sep=" ")
  #     }
  #     cat("\n")    
  #   }
  #   else if (is.data.frame(x$edges[[i]])) {
  #     for (j in seq_len(nrow(x$edges[[i]]))) {
  #       cat(x$vnames[x$edges[[i]][j,1]], edgeSymb[whEdge[i]],
  #                    x$vnames[x$edges[[i]][j,2]], "\n", sep=" ")
  #     }
  #     cat("\n")
  #   }      
  #   else if (is.list(x$edges[[i]])) {
  #     if (!is.null(x$edges[[i]]) && length(x$edges[[i]]) > 0) {
  #       for (j in seq_along(x$edges[[i]])) {
  #         cat(x$vnames[x$edges[[i]][[j]][1]], edgeSymb[whEdge[i]],
  #             x$vnames[x$edges[[i]][[j]][2]], "\n", sep=" ")
  #       }
  #       cat("\n")
  #     }
  #   }
  # } 
  # 
  invisible(x)
}

##' @export
`[.CADMG` <- function(graph, v, ..., drop=FALSE, etype) {
  ## 
  if (missing(v)) v <- graph$v
  else v <- unique.default(sort.int(v))
  
  if (!missing(etype)) gr_out <- as.mixedgraph(graph)[v, drop=drop, etype=etype]
  else gr_out <- as.mixedgraph(graph)[v, drop=drop]
  class(gr_out) <- c("CADMG", "mixedgraph")
  
  if (drop) gr_out$vtypes <- graph$vtypes[v]
  else {
    gr_out$vtypes <- graph$vtypes
    if (length(v) > 0) gr_out$vtypes[-v] <- NA
  }
  
  gr_out
}

##' @export
as.mixedgraph <- function(graph) {
  UseMethod("as.mixedgraph")
}

##' @export
as.mixedgraph.default <- function(graph) {
  convert(graph, format="mixedgraph")
}

##' @export
as.mixedgraph.CADMG <- function(graph) {
  class(graph) <- NULL
  graph <- graph[c("v","edges","vnames")]
  class(graph) <- "mixedgraph"
  graph
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
fix <- function(graph, w, v) {
  if (!is_ADMG(graph)) stop("Graph should be an ADMG")
  
  if (missing(v)) v = setdiff(graph$v, w)
  else if (length(intersect(w,v))) stop("Fixed and random vertex sets must be disjoint")
  else graph = graph[c(w,v)]
  
  v_left <- sort.int(c(w,v))
  
  if (any(!(v_left %in% graph$v))) stop(paste("No node ", paste(setdiff(c(w,v), graph$v), sep=", "), " found", sep=""))
  
  # graph <- withEdgeList(graph)
  
  graph <- mutilate(graph[v_left], w, etype="directed", dir=-1L)
  graph <- mutilate(graph, w, etype="bidirected", dir=0L)
  
  # e = graph$edges
  # if (nedge(graph, "directed") > 0) {
  #   
  #   gr_d <- graph[, etype="directed"]
  #   rmv = (matrix(unlist(e$directed), nrow=2)[2,] %in% w)
  #   e$directed = e$directed[!rmv]
  #   class(e$directed) = "eList"
  #   removeEdges(graph, )
  # }
  # if (nedge(graph, "bidirected") > 0) {
  #   bied = matrix(unlist(e$bidirected), nrow=2)
  #   rmvbi = (bied[1,] %in% w) | (bied[2,] %in% w)
  #   e$bidirected = e$bidirected[!rmvbi]
  #   class(e$bidirected) = "eList"
  # }
  # if (nedge(graph, "undirected") > 0) {
  #   stop("Only valid for CADMGs")
  # }

  vtype <- rep(NA_character_, length(graph$vnames))
  
  if (!is.null(graph$vtypes)) {
    vtype[v_left][graph$vtypes[v_left] == "fixed" | graph$v %in% w] <- "fixed"
    vtype[v_left][graph$vtypes[v_left] == "random" & graph$v %in% v] <- "random"
  }
  else {
    vtype[v_left][graph$v %in% v] <- "random"
    vtype[v_left][graph$v %in% w] <- "fixed"
  }
  
  graph$vtypes <- vtype
  
  out <- list(v=v_left, edges=graph$edges, vnames=graph$vnames, vtypes=vtype)
  class(out) <- c("CADMG", "mixedgraph")
  out
}

##' Get vector of fixed or random nodes in CADMG
##' 
##' @param graph a \code{mixedgraph} that is also a \code{CADMG}
##' 
##' @export
random <- function(graph) {
  if (is(graph, "CADMG")) {
    out <- na.omit(which(graph$vtypes == "random"))
    if (any(is.na(out))) stop("Failed to find random vertices")
    return(out)
  }
  else if (is.mixedgraph(graph)) return(graph$v)
  else stop("Not a valid object for random()")
}

##' @describeIn random get fixed nodes
##' @export
fixed <- function(graph) {
  if (is(graph, "CADMG")) {
    out <- na.omit(which(graph$vtypes == "fixed"))
    if (any(is.na(out))) stop("Failed to find fixed vertices")
    return(out)
  }
  else if (is.mixedgraph(graph)) return(integer(0))
  else stop("Not a valid object for fixed()")
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
