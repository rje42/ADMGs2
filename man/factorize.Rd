% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intrinsic.R
\name{partition}
\alias{partition}
\alias{factorize}
\title{Factorize a set of nodes into heads and tails}
\usage{
partition(graph, v = graph$v, r = TRUE, ht, head_order)

factorize(graph, v = graph$v, r = TRUE, ht, head_order)
}
\arguments{
\item{graph}{a CADMG}

\item{v}{integer vector of vertices to partition}

\item{r}{logical indicating whether nested parameterization is being used (used only if \code{ht} not provided)}

\item{ht}{the output of applying \code{headsTails()} to \code{graph}}

\item{head_order}{optional vector that represents partial order of heads}
}
\description{
Factorize a set of nodes into heads and tails
}
\section{Functions}{
\itemize{
\item \code{partition()}: Return the partition function for a particular set of vertices

}}
