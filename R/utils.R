## Authors 
## Benedikt Geier, bgeier@mail.uni-mannheim.de
## Borui Niklas Zhu, bzhu@mail.uni-mannheim.de
##
## Various help functions
##
## Copyright (C) 2019 -- 2020 Benedikt Geier
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


# given a family name (e.g. "beta"), a package name (e.g. "stats") 
# and the type ("r" for "rbeta", "d" for "dbeta") gets the function
# from the desired package and returns it
# this is needed as we cant use "stats::rbeta" directly in formals or 
# do.call -> error, that it cannot find the function

get_fun_from_package_internal <- function(type, fam, package) {
  return( get(paste0(type, fam), envir = asNamespace(package)) )
}

setGeneric(name = "get_fun_from_package",
           def = function(type, family, package) {
             standardGeneric("get_fun_from_package")
           })
setMethod(f = "get_fun_from_package",
          signature = c("character", "optimParams", "missing"),
          definition = function(type, family, package) {
            return(
              get_fun_from_package_internal(type, family@family, 
                                            family@package))
          })

setMethod(f = "get_fun_from_package",
          signature = c("character", "list", "missing"),
          definition = function(type, family, package) {
            return(get_fun_from_package_internal(type, family$family, 
                                                 family$package))
          })

setMethod(f = "get_fun_from_package",
          signature = c("character", "character", "character"),
          definition = function(type, family, package) {
            return(get_fun_from_package_internal(type, family, package))
          })

sample_params <- function(family, family_info, params = NULL, max = 100) {
  pars <- if (is.null(params)) family_info$lower else params
  
  # bound possible values to reasonable range
  upper <- pmin(family_info$upper, max)
  lower <- pmax(family_info$lower, -max)
  
  # simulate some "true" parameter values
  # cat("Sampling some 'true' parameter values and simulating sample data\n")

  float <- names(pars)[family_info$accepts_float]
  pars[float] <- runif(length(float), lower[float], upper[float])
  int <- names(pars)[!family_info$accepts_float]
  
  if (length(int) > 0) {
    impossible_ranges <- 
      (floor(family_info$upper[int]) - ceiling(family_info$lower[int]) < 0)
    if (any(impossible_ranges)) 
      stop("Param(s)",paste(names(int)[which(impossible_ranges)],collapse = ", "), 
           "do not accept floats andtheir range does not include any integer")
    pars[int] <- sapply(1:length(int), 
                        function(i) {
                          sample(ceiling(lower[int[i]]): floor(upper[int[i]]), 1)
                        }
		       )
  }

  # make sure that sampled values make sense for uniform distribution
  if (family$family == "unif") {
    h <- pars
    pars["min"] <- min(h)
    pars["max"] <- max(h)
  }
  return(pars)
}

sample_data <- function(n, family, params) {
  rfun <- get_fun_from_package(family = family$family, package = family$package, 
                               type = "r")
  n_or_nn <- if (! "nn" %in% names(formals(rfun))) list(n = n) else list(nn = n)
  args <- c(n_or_nn, params)
  testing_data <- do.call(rfun, args)
  return(testing_data)
}

# information criteria
informationCriteria <- function(ll, k, n) {
  aic = 2*k - 2 *  ll
  bic = log(n) * k - 2 *  ll
  aicc = aic + (2 * k^2 + 2 * k) / (n - k - 1)
  return(list(AIC = aic, BIC = bic, AICc = aicc))
}



is.natural <- function(x, tol = .Machine$double.eps^0.5) {
  (abs(Im(x)) < tol) && (abs(Re(x)) > tol) &&
  isTRUE(all.equal(x, round(x),
                   tolerance = tol,
                   check.attributes = FALSE,
                   check.names = FALSE))
}

install.packages_DistributionFitR <- function(...) {
  if(hasArg("pkgs")) stop("Do not specify argument 'pgks'")
  installed <- rownames(installed.packages())
  all <- unique(sapply(FamilyList, function(x) x$package))
  needed <- setdiff(all, installed)
  install.packages(needed, ...)
}
