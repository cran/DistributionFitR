\name{globalfit}
\alias{globalfit}

\title{
  Detect continuity and fit multiple distributions to given data
}
\description{
  Given a numerical data vector, this function fits multiple distributions 
  with the maximum likelihood method and returns an object containing the best 
  fitted parameters and information criteria.
  Refer to the \dQuote{Examples} section or the result class
  \code{\link[=globalfit-class]{globalfit}} on how to sort and output the results with
  e.g. \code{\link[=summary,globalfit-method]{summary}}.
}
\usage{
globalfit(data, continuity = NULL, method = "MLE",
  verbose = TRUE, packages = "stats",
  append_packages = FALSE, cores = NULL,
  max_dim_discrete = Inf, sanity = 1,
  timeout = 5
  )
}

\arguments{
  \item{data}{
    numeric vector of data points.
  }
  \item{continuity}{
    logical; if \code{TRUE}, the data is fitted with continuous distributions.
    If no input is given, the data will be tested for continuity.
  }
  \item{method}{
    character;
    method for parameter estimation. So far only Maximum-Likelihood is 
    implemented, thus this argument must
    be \code{"MLE"}.
  }
  \item{verbose}{
    logical;
    if \code{TRUE}, show progress and packages from where to fit distributions.
  }
  \item{packages}{
     either a character vector with names of packages; \cr
    or a list such as those returned by \command{getFamilies}\cr
    or \code{NULL}, i.e. all families known by this package (recommended).  

    default: \dQuote{stats}. \cr
 } 
  \item{append_packages}{
   logical; if \code{TRUE} (default) appends packages specified in the 
   argument \code{packages} to the standard search list, if \code{FALSE} 
   \code{globalfit} will use only those packages and ignore the standard
   search list.
  }
  \item{max_dim_discrete}{
   non-negative integer; distributions with more 
   non-continuous parameters than \code{max_dim_discrete} will not be 
   considered. Manual setting is recommended if calculation speed has to
   be cut down.
  }
  \item{cores}{
   integer; number of CPU cores to be used in the calculations of best fitted 
   parameters and information criteria.
  }
  \item{sanity}{
    either a positive numeric or logical;
    if it is a positive numeric, it
    controls a sanity check where obviously bad fits
   are filtered out. The smaller the number, the stricter the check will be
   executed and the more potential distributions will be rejected.

   If \code{sanity = FALSE} a sanity check is not carried out. \cr
   (DistributionFitR generally depends on other packages to supply reasonable 
   distribution functions.)

   Default is 1. 
 }
 \item{timeout}{
   logical or numeric.
   if it is a
   positive numeric, it gives the
   seconds until timeout for the underlying optimiser
    \code{\link[stats]{optim}}.

    If \code{timeout = FALSE} no timeout is performed.
   }
}

\details{
  If there is no continuity input given, this function first tests via multiple 
  criteria whether the data is continuously or discretely distributed. 
  Given that information, the related distributions from \code{getFamilies()} 
  are fitted to the data via maximum likelihood method and information 
  criteria are calculated. For discrete data not in the form of 
  integers only, an appropriate linear transformation is applied to ensure
  stable optimization.

  Since DistributionFitR technically allows for comparing over all distributions
  in any R-package, computation speed is likely to be an issue. The
  following may 
  help:
  \itemize{
    \item{using argument \code{packages} with \code{append_packages = FALSE}
	    to restrict the search to certain packages
    }
    \item{discarding distributions with too many discrete parameters using
	    argument \code{max_dim_discrete}
    }
    \item{specifying \code{timeout}, which affects the maximum time spent on
	    each distribution (not overall!). The value in \code{timeout}
	    will not be translated directly to the actual maximum time
	    due to differing number of times  \code{optim} is run
	    under different algorithms.
    }
  }
}

\value{
  \code{globalfit} returns an object of class
  \code{\link[=globalfit-class]{globalfit}}.
}

\author{
  Moritz Lauff, Kiril Dik, Nadine Tampe, Borui Niklas Zhu, 
  Benedikt Geier, Moritz Kern
}

\examples{
  # Example 1
  data <- rnorm(n = 100, mean = 70, sd = 4)
  r <- globalfit(data, cores = if(interactive()) NULL else 2)
  summary(r)

  # continuous or discrete
  \donttest{
  # Example 2
  # Alternatively, it is possible to input whether the data is
  globalfit(data, continuity = TRUE)

  # Example 3
  # fit over all distribution in the standard search list
  globalfit(data, packages = NULL)
  }
}
\keyword{distribution}
