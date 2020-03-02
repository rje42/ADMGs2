#' Summarize fitted ADMG model.
#' 
#' Produces summary statistics for fit of data to an ADMG model.
#' 
#' @aliases summary.mixed_fit print.mixed_fit_summary
#' @param object An object of class \code{mixed_fit}.
#' @param fisher logical: should Fisher Information Matrix be calculated and
#' standard errors provided?
#' @param \dots Other arguments to pass to \code{summary}.
#' @return List containing the output from \code{\link{fitADMG}} as well as:
#' \item{p}{The number of parameters.} \item{AIC}{Akaike Information Criterion
#' of model fit.} \item{BIC}{Bayesian Information Criterion of model fit.}
#' \item{deviance}{The deviance of the model fit.} \item{n}{The dimension of
#' the data.} \item{se.table}{A table containing standard errors (if
#' \code{fisher = TRUE}).} \item{FIM}{The Fisher Information Matrix (if
#' \code{fisher = TRUE}).} \item{probs}{Best fitting probability distribution.}
#' @author Robin Evans
#' @seealso \code{\link{fitADMG}}.
#' @references Evans, R.J. and Richardson, T.S. (2010) - Fitting acyclic
#' directed mixed graphs to binary data. \emph{UAI-10}.
#' @keywords graphs
#' 
#' @method summary mixed_fit
summary.mixed_fit <-
function (object, fisher=TRUE, ...) {
  q = object$params$q
  p = length(unlist(q)) # i.e. the numbers of parameters
  n_obs = sum(object$dat)
  
  AIC = 2*p - 2*object$ll
  BIC = log(n_obs)*p - 2*object$ll
  
  if (is.array(object$dat)) ll.sat = sum(object$dat*log(object$dat/n_obs))
  else ll.sat = sum(object$dat$freq*log(object$dat$freq/n_obs))
  deviance = 2*(ll.sat - object$ll)
  
  if (fisher) FIM = .fisher_mixed_fit(object)
  else FIM = NULL
  probs = probdist(object$params, object$maps)

  se.table = matrix(NA, ncol=2, nrow=p, dimnames = list(getMparamsNames(object$params), c("Estimate", "Std. Error")))
  qvec = unlist(q)
  se.table[,1] = qvec
  if (!is.null(object$SEs)) se.table[,2] = object$SEs
  else if (fisher) se.table[,2] = sqrt(diag(solve.default(FIM)))
  else se.table[,2] <- NA

  out = c(object, list(p=p, AIC=AIC, BIC=BIC, deviance=deviance, se.table=se.table, FIM=FIM, probs=probs))
  class(out) = "mixed_fit_summary"
  out
}


##' @method print mixed_fit
print.mixed_fit <-
  function (x, ...) {
    q = x$params$q
    p = length(unlist(q)) #- length(q)
    
    cat("Acyclic Directed Mixed Graph fit to discrete data\n")
    cat("Variable dimensions:", x$dim, "\n", sep=" ")
    cat(ifelse(x$r, "Recursive", "Non-recursive"), "parametrization\n", sep=" ")
    cat("log-likelihood: ", x$ll, " on ", p, " parameters\n")
    invisible(x)
  }


##' @method print mixed_fit_summary
print.mixed_fit_summary <-
  function (x, ...) {
    cat("Acyclic Directed Mixed Graph fit to discrete data\n")
    cat("Variable dimensions:", x$dim, "\n", sep=" ")
    cat(ifelse(x$r, "Recursive", "Non-recursive"), "parametrization\n\n", sep=" ")
    
    cat("log-likelihood: ", x$ll, " on ", x$p, " parameters\n\n", sep="")
    
    cat("AIC = ", x$AIC, "\n", sep="")
    cat("BIC = ", x$BIC, "\n", sep="")
    cat("Deviance = ", x$deviance, " on ", prod(x$dim)-x$p-1, " degrees of freedom\n", sep="")
    
    cat("\n")
    print(signif(x$se.table,4))
    invisible(x)
  }

