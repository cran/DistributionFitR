\name{optimParams-class}
\alias{optimParams-class}

\title{Class \dQuote{optimParams}}
\description{The class \code{optimParams} handles the slot
	\code{fits} of the class \code{\link[=globalfit-class]{globalfit}}.
	It contains for some given data the optimization results
	and their characteristics.}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("optimParams", family, package, estimatedValues, AIC, BIC, AICc, 
             sanity)}.
  More comfortably, you may use the function \code{\link{globalfit}}.
  The list elements in the \code{fits}-slot of the result is a
  \code{optimParams}-object.
}

\section{Slots}{
  \describe{
    \item{family}{character string; indicating the family name in the package,
                  e.g. \dQuote{norm}.}
    \item{package}{character string; indicating the package where the family 
                   was found, e.g. \dQuote{stats}.}
    \item{estimatedValues}{named numeric vector; the estimated 
                           parameters of the family.}
    \item{AIC}{numeric; Akaike Information Criterion.}
    \item{BIC}{numeric; Bayes Information Criterion.}
    \item{AICc}{numeric; small sample-corrected Akaike Information Criterion.}
    \item{sanity}{named list; with numeric items \code{hist_check} und 
                  \code{int_check} and boolean item \code{good} as a 
                  result of the sanity check.}
  }}
  

\author{
  Moritz Kern
}

\seealso{
  \code{\link[=globalfitSummary-class]{globalfitSummary}},
  \code{\link{globalfit}}
}

\examples{
  data <- rnorm(n=100, mean=10, sd= 1)
  r <- globalfit(data, cores= if (interactive()) NULL else 1,
                 packages="stats", append_packages=FALSE)
  
  # optimParams of best fit
  r@fits[[1]]
    
}

\keyword{distribution}
