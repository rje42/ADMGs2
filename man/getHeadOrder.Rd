% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maps.R
\name{getHeadOrder}
\alias{getHeadOrder}
\title{Get head order from graph}
\usage{
getHeadOrder(graph, heads, check = FALSE, r = FALSE)
}
\arguments{
\item{graph}{an object of class \code{mixedgraph}}

\item{heads}{list of heads}

\item{check}{logical: should we check that each element of \code{head} could actually be a head?}

\item{r}{logical: should heads be treated as recursive?}
}
\description{
Get head order from graph
}
\details{
Returns integer vector giving a partition suitable 
ordering for heads.  Not certain this works if \code{r=TRUE}.
}
