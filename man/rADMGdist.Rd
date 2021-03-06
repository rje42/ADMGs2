% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probdist.R
\name{rADMGdist}
\alias{rADMGdist}
\title{Generate a distribution from an ADMG model}
\usage{
rADMGdist(graph, dims, map, r = TRUE, alpha = 1)
}
\arguments{
\item{graph}{an object of class graph}

\item{dims}{integer vector of dimensions of distribution}

\item{map}{optionally, the output of \code{maps()}}

\item{r}{logical indicating whether or not recursive factorizations should
be used}

\item{alpha}{parameter to scale Beta distribution down by.}
}
\value{
Object of class \code{moebius} giving generalized Moebius parameters
for a distribution in the model associated with \code{graph}.
}
\description{
Generate a distribution from an ADMG model
}
\details{
The random distribution is
obtained by starting with the uniform distribution, and adding an 
independent (scaled) normal random variables to each Moebius parameter.  We 
then scale the move to reach approximately the boundary of the simplex,
and then pick a value between the uniform distribution and the boundary 
point using a Beta distribution.  The parameter \code{alpha} being larger
favours distributions closer to the uniform distribution.
}
