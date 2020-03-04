\name{getFamily}
\alias{getFamily}

\title{
   Find all distributions in a package
}
\description{
   Given the name of a package, this function finds all distribution families
   that are provided in the package.
   Distributions are identified by scanning all function names in a package
   for the pattern r***, p***, q***, d***,
   (*** representing at least one character)
}
\usage{
getFamily(pkg)
}

\arguments{
  \item{pkg}{
    character string; name of the package
  }
}

\note{
  The function make do with
  at least two functions that have the same remainder of the function
  name.
}

\value{
  A list of of lists, the latter with two named elements:
  \item{package}{The name of the package as provided to the function}
  \item{family}{The name of the function that belongs to a distribution family,
  referred to as *** above.}
}
\author{
   Manuel Hentschel, Valentin von Trotha
}

\examples{
   str(getFamily("stats"))
}

\keyword{distribution}
