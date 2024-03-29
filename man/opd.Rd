% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opd.R
\name{opd}
\alias{opd}
\title{Get the ordered power DAG}
\usage{
opd(graph, topOrd, v, include_empty = FALSE)
}
\arguments{
\item{graph}{an ADMG of class \code{mixedgraph}}

\item{topOrd}{a topological ordering of the vertices in \code{graph}}

\item{v}{optionally the vertices to obtain components for}

\item{include_empty}{logical: should independences be included even if there are no variables in one set?}
}
\description{
Returns the ordered power DAG for an ordinary Markov ADMG
}
\details{
Returns the ordered power DAG, as defined in Hu and Evans (2022),
giving independences implied by the local Markov property and sufficient to
define the model.
}
\references{
Z. Hu and R.J. Evans. Towards standard imsets for maximal ancestral
graphs, _arXiv:2208.10436_, 2022.
}
