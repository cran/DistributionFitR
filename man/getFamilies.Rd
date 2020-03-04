\name{getFamilies}
\alias{getFamilies}
\title{
   Lists (all) distribution families with their parameters
}
\description{
   This function extracts distribution families from R-packages along
   with their  
   parameters and characteristics, and scans specific packages if desired.
}
\usage{
getFamilies(all.packages)
}

\arguments{
  \item{all.packages}{
    character vector or missing or logical; package names in which to scan:
    \describe{
      \item{character}{
	 the distribution families given
      }
      \item{\code{TRUE}}{all installed packages}
      \item{\code{FALSE}}{base packages of R}
    }
    If missing, a list of recognized families is returned.
  }
}

\details{
  This function retrieves distribution families from R-packages as specified 
  in \sQuote{arguments}. A distribution family is defined as follows: \cr
  Take the functions \code{dnorm}, \code{pnorm}, \code{qnorm} and \code{rnorm} 
  from the package \sQuote{stats}. These functions are related to sampling, 
  density etc. of the normal distribution. The family name is the part of 
  the function name without the prefixes \sQuote{d},  \sQuote{p}, \sQuote{q}, or
  \sQuote{r} , in this case: \sQuote{norm}.
}
\value{
  A list of lists. Each sublist pertains to exactly one distribution 
  family (such as \sQuote{norm} from \sQuote{stats}) and contains:
  \item{package}{
    character; name of the package containing the family
  }
  \item{family}{
    character; name of the family, as defined under \dQuote{Details}
  }
  \item{family_info}{
    \itemize{
      \item{\code{lower}:}{
	named numeric vector; lower bounds for distribution parameters
      }
      \item{\code{upper}:}{
	named numeric vector; upper bounds for distribution parameters
      }
      \item{\code{accepts_float}:}{
	named logical vector; whether each parameter can be any real 
	number (\code{TRUE}) or only integers are valid inputs (\code{FALSE})
      }
      \item{\code{defaults}:}{
	names numeric vector; set of default parameters that jointly constitute valid 
	input. Used for optimisation.
      }
      \item{\code{log}:}{
	logical; whether log values are provided by the density function
      }
      \item{\code{discrete}:}{
	logical; whether the distribution itself (not the parameters) takes on 
	discrete values only
      }
      \item{\code{support_min}:}{
	numeric; lower bound of the support, i.e. lowest value of \eqn{x} 
	where \eqn{f(x) > 0} numerically for any valid parameter set, \eqn{f} the density function
      }
      \item{\code{support_max}:}{
	numeric; upper bound of the support
      }
      \item{\code{support_max_depends_on}:}{
	named logical vector; whether \code{support_max} depends on the value of each 
	parameter
      }
      \item{\code{support_min_depends_on}:}{
	named logical vector; whether \code{support_min} depends on the value of each 
	parameter
      }
    }
  }
}

\author{
   Tim Glockner, Adrian Heppeler, Borui Niklas Zhu
 }

\examples{ 
  str(getFamilies())
}
\keyword{distribution}
