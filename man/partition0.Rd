% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intrinsic.R
\name{partition0}
\alias{partition0}
\alias{factorize0}
\title{Gives partition/factorization}
\usage{
partition0(graph, heads, v = seq_len(graph$n), r = TRUE, head.order)

factorize0(graph, v = seq_len(n), r = TRUE, ht, head.order)
}
\arguments{
\item{graph}{object of class \code{mixedgraph}}

\item{heads}{list of heads}

\item{v}{set of vertices to partition}

\item{r}{logical: should recursive parameterization be used?}

\item{head.order}{numeric vector in same order as heads}

\item{ht}{list of heads and tails}
}
\description{
Uses order for speed
returns integer value of heads from provided list
}
\section{Functions}{
\itemize{
\item \code{factorize0()}: Give factorization of heads and tails

}}
