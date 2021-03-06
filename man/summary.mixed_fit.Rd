% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixed_fit.R
\name{summary.mixed_fit}
\alias{summary.mixed_fit}
\alias{print.mixed_fit_summary}
\title{Summarize fitted ADMG model.}
\usage{
\method{summary}{mixed_fit}(object, fisher = TRUE, ...)
}
\arguments{
\item{object}{An object of class \code{mixed_fit}.}

\item{fisher}{logical: should Fisher Information Matrix be calculated and
standard errors provided?}

\item{\dots}{Other arguments to pass to \code{summary}.}
}
\value{
List containing the output from \code{\link{fitADMG}} as well as:
\item{p}{The number of parameters.} \item{AIC}{Akaike Information Criterion
of model fit.} \item{BIC}{Bayesian Information Criterion of model fit.}
\item{deviance}{The deviance of the model fit.} \item{n}{The dimension of
the data.} \item{se.table}{A table containing standard errors (if
\code{fisher = TRUE}).} \item{FIM}{The Fisher Information Matrix (if
\code{fisher = TRUE}).} \item{probs}{Best fitting probability distribution.}
}
\description{
Produces summary statistics for fit of data to an ADMG model.
}
\references{
Evans, R.J. and Richardson, T.S. (2010) - Fitting acyclic
directed mixed graphs to binary data. \emph{UAI-10}.
}
\seealso{
\code{\link{fitADMG}}.
}
\author{
Robin Evans
}
\keyword{graphs}
