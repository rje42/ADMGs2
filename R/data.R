##' Data Sets
##' 
##' Some datasets for illustration.
##' 
##' @format These are all \code{data.frame} objects.
##' 
"twins"

##' @describeIn twins Data from the General Social Survey
"gss"

##' @describeIn twins Subset of data from the General Social Survey
"gss_small"


##' van Alstyne and Gottfredson Data
##'
##' A criminology dataset for illustration.
##'
##' @format This is a \code{data.frame} object, with the following
##' variables:
##' \describe{
##'   \item{person}{Was this a crime against a person?}
##'   \item{age}{Was the perpatrator aged over 25?}
##'   \item{drug}{Did they have a drug of alcohol dependancy?}
##'   \item{validation}{Indictator that this was from a validation sample}
##'   \item{success}{Was the parolee successful in avoiding arrest for further offences?}
##'   \item{prior}{Was there a previous sentence for this offence?}
##'   \item{freq}{The total number of individuals in this row}
##' }
##'
##' @source van Alstyne and Gottfredson, A Multidimensional Contingency Table Analysis of Parole Outcome: New Methods and Old Problems in Criminological Prediction,
##' \emph{J. Res. Crime and Deliquency}, 1978.
"alstyne"


##' Caesarian Dataset
##'
##' A dataset on Caesarian sections for illustration.
##'
##' @format This is a \code{data.frame} object, with the following
##' variables:
##' \describe{
##'   \item{cesarian}{Was a Caesarian section performed?}
##'   \item{monitor}{}
##'   \item{arrest}{}
##'   \item{breech}{Was the baby in the breech position?}
##'   \item{nullipar}{Is this a first birth?}
##'   \item{freq}{The total number of births in this row}
##' }
##'
##' @source Unknown.
"caesarian"

##' Sewell and Shah Dataset
##'
##' A dataset on educational aspirations
##'
##' @format This is a \code{data.frame} object, with the following
##' variables:
##' \describe{
##'   \item{Sex}{0=female, 1=male}
##'   \item{Iq}{IQ, rated from 0 to 3}
##'   \item{Cp}{Indicator of college plans}
##'   \item{Pe}{Indicator of parental encouragement}
##'   \item{Ses}{Socio-economic status, from 0 to 3}
##' }
##'
##' @source Sewell, W. H., & Shah, V. P. (1968). Social class, parental encouragement, and educational aspirations. \emph{American Journal of Sociology}, 73(5), 559-572.
##' @seealso \url{https://doi.org/10.1086/224530}
"educ_asp"

##' The '16-not-18' Graph
##' 
##' The '16-not-18 graph' is an example of a graph with a nested constraint
##' the removes two parameters in the all binary case: hence the name 16 
##' (parameters) not 18 (parameters).  
##' 
##' @format A mixed graph object that is also an ADMG.
##' 
"gr16n18"

##' Graph that often serves as a counter example
##' 
##' @format A mixed graph object that is also an ADMG.
##' 
"grCounter"