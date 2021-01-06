##' @export
districts.CADMG <- function(graph, ...) {
  gr0 <- graph[random(graph)]
  gr0 <- mixedgraph(v=random(graph), edges = gr0[etype="bidirected"]$edges, 
             vnames=vnames(graph))
  districts(gr0)
}