\name{globalfitSummary-class}
\alias{globalfitSummary-class}
\alias{show,globalfitSummary-method}
\alias{print,globalfitSummary-method}

\title{Class "globalfitSummary"}
\description{The \code{globalfitSummary}-class is the class to handle the 
  entries of the class
  \code{\link[=globalfit-class]{globalfit}}
  and assemble them together
  for summary.}
  
\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("globalfitSummary", data, continuity, method, fits, ic)}.
  More comfortably, you may use the method
  \code{\link[base]{summary}}
  on an object of class \code{\link[=globalfit-class]{globalfit}} -
  its result is a \code{globalfitSummary}-object.
}

\section{Slots}{
  \describe{
    \item{call}{the call, which created the 
    \code{\link[=globalfit-class]{globalfit}}, where this object originated.}
    \item{data}{vector of data points}
    \item{continuity}{logical; if \code{TRUE}, indicates that the data 
          points come from a continuous distribution;
	  if \code{FALSE}, indicates that they come from a
	  discrete distribution}
    \item{method}{character; the method used for the fit.}
    \item{fits}{data frame; sorted by the \code{ic} selected in
	    \code{summary} or the constructor call, 
	    with the columns: family, package and ic.}
    \item{ic}{character; indicates by which crietrion \code{fits} was sorted.}    
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(x = "globalfitSummary")}: display the object}
    \item{print}{\code{signature(x = "globalfitSummary")}: calls \code{show}}
    }
}
  
\author{
  Moritz Kern
}

\seealso{
  \code{\link[=globalfitSummary-class]{globalfitSummary}} \cr
  \code{\link{globalfit}} \cr
  \code{\link[=optimParams-class]{optimParams}}
}

\examples{
  data <- rnorm(n=100, mean=10, sd= 1)
  r <- globalfit(data, cores= if (interactive()) NULL else 1,
                 packages="stats", append_packages=FALSE)
  summary(r)
}

\keyword{distribution}
