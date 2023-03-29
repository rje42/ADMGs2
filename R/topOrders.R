##' @export
all_tops <- function (graph) {
  PAs <- lapply(graph$v, function(x) pa(graph,x))
  CHs <- lapply(graph$v, function(x) ch(graph,x))
  
  inDeg <- lengths(PAs)
  sq <- seq_along(graph$v)
  
  out <- list()
  vst <- rep(FALSE, length(CHs))
  
  ## get vertices with no parents
  zeroIn <- which(inDeg == 0)
  if (length(zeroIn) == 0) stop("Graph is cyclic")
  
  for (k in seq_along(zeroIn)) {
    inDeg[sq[CHs[[sq[k]]]]] <- inDeg[sq[CHs[[sq[k]]]]] - 1
    vst[sq[k]] <- TRUE
    
    tmp <- to_recurse(CHs, inDeg, vst)
    tmp <- lapply(tmp, function(x) c(graph$v[k], x))
    
    inDeg[sq[CHs[[sq[k]]]]] <- inDeg[sq[CHs[[sq[k]]]]] + 1
    vst[sq[k]] <- FALSE
    
    out <- c(out, tmp)
  }

  out
}

to_recurse <- function (CHs, inDeg, vst) {
  
  if (all(vst)) return(list(integer(0)))

  ## get vertices with no parents
  zeroIn <- which(inDeg == 0 & !vst)
  if (length(zeroIn) == 0) stop("Graph is cyclic")
  
  out <- list()
  
  for (k in zeroIn) {
    inDeg[CHs[[k]]] <- inDeg[CHs[[k]]] - 1
    vst[k] <- TRUE
    
    tmp <- to_recurse(CHs, inDeg, vst)
    tmp <- lapply(tmp, function(x) c(k, x))
    
    inDeg[CHs[[k]]] <- inDeg[CHs[[k]]] + 1
    vst[k] <- FALSE
    
    out <- c(out, tmp)
  }
  
  out  
}
