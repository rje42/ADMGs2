% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsets.R
\name{subsetRep}
\alias{subsetRep}
\title{Get the subset representation of an ordinary ADMG}
\usage{
subsetRep(graph, max_size, sort = 1, r = FALSE)
}
\arguments{
\item{graph}{ADMG}

\item{max_size}{largest set to consider}

\item{sort}{if 1, returns unique sets, if 2 returns sorted sets, if 3 returns sets ordered reverse lexicographically}

\item{r}{logical: use recursive heads? (Defaults to \code{FALSE})}
}
\description{
Get the subset representation of an ordinary ADMG
}
\details{
This returns sets of the form \code{HuA}, where \code{A} is any 
subset of \code{tail(H)}.  This representation is shown to 
correspond precisely to ordinary Markov equivalence.
}
