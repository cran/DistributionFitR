\name{install.packages_DistributionFitR}
\alias{install.packages_DistributionFitR}

\title{
  Installs all packages from DistributionFitR's standard search list.
}
\description{
  DistributionFitR comes with an extensive list of distribution families
  on CRAN together with their characteristics, referred to in our manuals
  as \dQuote{standard search list}. This function is a helper to install all
  the packages which contribute to this list. Usage is exactly like
  \code{install.packages} in base R, but without the need to specify
  which packages to install. If some packages fail to install, the others
  will not be affected (except dependencies).
}
\usage{
  install.packages_DistributionFitR(...)
}

\arguments{
  \item{...}{
    any argument that can be passed to base R's \code{install.packages}, except
    the argument \code{pkgs}.
  }
}

\value{
  invisible \code{NULL}.
}
\author{
  Borui Niklas Zhu
}

\examples{
  \dontrun{
   # running the next command will run several minutes
   # and install many packages
   install.packages_DistributionFitR()
  }
}

\keyword{install}
\keyword{packages}
