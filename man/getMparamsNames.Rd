% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maps.R
\name{getMparamsNames}
\alias{getMparamsNames}
\title{Return Moebius parameter name strings}
\usage{
getMparamsNames(mparam, blanks = FALSE)
}
\arguments{
\item{mparam}{An object of class \code{mparam}.}

\item{blanks}{logical indicating whether to omit repeated names.}
}
\value{
A string valued vector containing the names in head order.
}
\description{
Gives a vector of strings corresponding to the meaning of each generalized
Moebius parameter: e.g. "p(a, c = 0 | b = 1, d = 0, e = 0)".
}
\details{
The logical \code{blanks} provides for adding whitespace in place of some
repeated vertices in the output.
}
\note{
Note that if the Moebius parametrization is with respect to the
recursive factorization, some of the parameters may be reweighted
probabilities of the joint distribution.
}
\examples{

data(gr2, package="MixedGraphs")
mparams = moebius(gr2, prob_table(gr2, values=rep(1/32, 32)))

ADMGs2:::getMparamsNames(mparams)


}
\references{
Evans, R.J. and Richardson, T.S. (2010) - Fitting acyclic
directed mixed graphs to binary data. \emph{UAI-10}.
}
\author{
Robin Evans
}
