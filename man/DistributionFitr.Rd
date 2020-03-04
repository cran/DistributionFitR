\name{DistributionFitR-package}
\alias{DistributionFitR}
\alias{DistributionFitR-package}
\docType{package}
\title{
  Fitting Multiple Parametric Distributions
}
\description{
Parametric densities (or count densities) will be fitted to user-given
univariate data via maximum likelihood. The user more or less only enters the
data. The program automatically searches for parametric distributions and
parameters thereof that best describe the data. It then returns 
the ten best parameter families including the fitted parameters. 

DistributionFitR comes in with a standard search list of 408 parametric
distribution families as given by R-packages on CRAN.

The package contains the following functions:
\itemize{
  \item \code{\link{globalfit}} Given some univariate data, determines the best 
  fitting parametric distributions from R-packages
  \item \code{\link{install.packages_DistributionFitR}} Installs R-packages
  that contribute to DistributionFitR's standard search list
  \item \code{\link{getFamilies}} Browsing for distribution families
  \item \code{\link{getFamily}} Find all distributions within a package
  \item \code{\link{getParams}} For a single distribution family, parameters are 
determined; upper/lower limits and other chraracteristics are extracted
%\item \code{optimParamsDiscrete}: Fits the parameters for a given
%distribution family with known names for the parameters and known value
%range of each parameter item 
%\item \code{loglik} Returns the loglikelihood function to be optimized item
}
}

\note{
The most relevant function is \code{\link{globalfit}}, with an S4 object
as return value. See the examples on how to display the results, such as
using \code{\link[base]{summary}} or
\code{\link[graphics]{hist}}.
Read more on the return value in
\code{\link[=globalfit-class]{globalfit}} and
\code{\link[=globalfitSummary-class]{globalfitSummary}}.

For exotic packages used frequently or where parameter extraction is
time-consuming, users may be interested to do the latter once with
\code{\link{getFamilies}} and save the results for subsequent usage as
argument in \code{\link{globalfit}}. The functions
\code{\link{getFamily}} and \code{\link{getParams}} are lower-level
functions invoked by \code{\link{getFamilies}}, and may be of usage for
other purposes. 
}

\author{
Borui Niklas Zhu, Benedikt Geier, Moritz Kern, Kiril Dik, Moritz Lauff, 
Manuel J. Hentschel, Adrian Heppeler, Niclas Lietzow, Till Freihaut, 
Tim Glockner, Nadine Tampe, Leonardo Vela, Helene Peter, 
Martin Schlather, Yiqi Li
  
Maintainer: Borui Niklas Zhu, \email{bzhu@mail.uni-mannheim.de}   
}

\seealso{
 \code{\link{globalfit}},
 \code{\link{install.packages_DistributionFitR}},
 \code{\link{getFamilies}},
 \code{\link{getFamily}},
% \code{\link{loglik}},
% \code{\link{optimParamsDiscrete}},
 \code{\link{getParams}}
}

\examples{ 
# example for globalfit
data <- rnorm(n = 100, mean = 70, sd = 4)
r <- globalfit(data, cores = if(interactive()) NULL else 2)
summary(r)

# example for getFamily
str(getFamily("stats"))

# example for getParams
getParams("beta", package = "stats")

}


\keyword{distribution}
\keyword{parameters}
\references{
  \itemize{
    \item \litZucchini
    \item \litShao
  }
}
