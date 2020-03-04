\name{globalfit-class} 
\alias{globalfit-class} 
\alias{sort,globalfit-method}
\alias{summary,globalfit-method}
\alias{print,globalfit-method}
\alias{AIC,globalfit-method}
\alias{BIC,globalfit-method}
\alias{hist,globalfit-method}

\title{Class \dQuote{globalfit}}
\description{The class \code{globalfit} handles return objects 
  from \link{globalfit}. It contains for some given data a list of fitted 
  distributions, their estimated parameters and supplementary information.}
\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("globalfit", data, continuity, method, fits)}.
  More comfortably, you may use the  function \code{\link{globalfit}}.
  The result of these calls is a globalfit object.
}

\section{Slots}{
  \describe{
    \item{call}{the call, which created this object}
    \item{data}{vector of data points}
    \item{continuity}{logical; if \code{TRUE}, indicating that the data 
          points come from a continuous distribution; 
	  if \code{FALSE}, indicating that they come from a discrete
	  distribution}
    \item{method}{character; the method used for the fit.}
    \item{fits}{list of S4-objects of class
      \code{\link[=optimParams-class]{optimParams}} 
  }}
}

\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "globalfit")}: 
                        summarizes the object and creates an object of 
                        \code{\link[=globalfitSummary-class]{globalfitSummary}}.
			Specify argument \code{ic} to choose how the results
			are to be sorted (as in method \code{sort}.

                  }
    \item{hist}{\code{signature(x = "globalfit")}: 
                computes a histogram of the given data points and plots it
		together with the density of the estimated best fit.
		Specify argument \code{which} to choose which fitted density
    		to overlay: the number of the fit as returned by
    		\code{summary}; i.e. \code{which = 1} for the best fit,
		\code{which = 2} for the second-best etc. \cr
    		Default is \code{1}.}
    \item{print}{\code{signature(x = "globalfit")}: applies the
 method \code{\link[base]{summary}} and prints the result.}
%    \item{IC}{\code{signature(x = "globalfit")}: returns the selected 
%                                                 information
%              criteria of the items in \code{fits} slot}
    \item{AIC}{\code{signature(x = "globalfit")}: shows the AIC value of the
	       fits. Specify argument \code{n} to display AIC for the \code{n}
	       best fits according to this criterion.}
    \item{BIC}{\code{signature(x = "globalfit")}: shows the BIC value of the
	       fits. Specify argument \code{n} to display BIC for the \code{n}
	       best fits according to this criterion.}
%    \item{AICc}{\code{signature(x = "globalfit")}: applies the \code{IC}-method 
%               for the small sample-adjusted Akaike Information Criterium.}
    \item{sort}{\code{signature(x = "globalfit")}: sorts the results in slot 
                \code{fits} by the information criterium selected.
		in argument \code{ic}.
                Available options are \code{"AIC"}, \code{"BIC"} or
		\code{"AICc"}.
    	       }
  }
}
  
\author{
  Moritz Lauff, Kiril Dik, Nadine Tampe, Borui Niklas Zhu, 
  Benedikt Geier, Moritz Kern
}

\seealso{
  \code{\link[=globalfitSummary-class]{globalfitSummary}} \cr
  \code{\link[=optimParams-class]{optimParams}} \cr
  \code{\link{globalfit}}
}

\examples{
  data <- rnorm(n = 100, mean = 10, sd = 1)
  r <- globalfit(data, cores = if(interactive()) NULL else 2)
                 
  sort(r, ic = 'BIC')
    
  print(r)
  
  summary(r)
  summary(r, ic = 'AICc', n = 7)
      
  hist(r, ic = 'BIC', which = 4)
  
  AIC(r, n = 2)
  BIC(r)
}

\keyword{distribution}
