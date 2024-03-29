% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probdist.R
\name{rDAGdist}
\alias{rDAGdist}
\title{Generate a distribution from a DAG model}
\usage{
rDAGdist(n, graph, dims = rep(2L, nv(graph)), alpha = 1)
}
\arguments{
\item{n}{number of distributions to simulate}

\item{graph}{a DAG object of class \code{mixedgraph}}

\item{dims}{integer vector of dimensions of distribution}

\item{alpha}{parameter to use in Dirichlet distribution}
}
\description{
Simulate conditional probability tables to obtain a 
discrete DAG model that is Markov with respect to the
graph.
}
\details{
Simulated conditional probability tables for each
node given its parents, and uses this to obtain the joint 
distribution.
}
