## Authors 
## Moritz Kern, mkern@mail.uni-mannheim.de
##
## output methods for globalfit-objects
##
## Copyright (C) 2019 -- 2020 Moritz Kern
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


"%@%" <- function(x, ic) eval(parse(text = paste0('x@', ic)))
  
setMethod(f = "sort", signature = c('globalfit'),
          def = function(x, decreasing = FALSE, ic = c('BIC', 'AIC', 'AICc')) {
            ic <- match.arg(ic)
            ic <- sapply(x@fits, function(x) x %@% ic) 
            x@fits <- x@fits[order(ic)]
            return(x)
          })



setMethod(f = "summary", signature = c("globalfit"),
          def = function(object, n = 10,  ic = c('BIC', 'AIC', 'AICc')) {
            ic <- match.arg(ic)
            
      	    if(length(object@fits) < 1) {
      	      warning("No family fitted. Either packages provided in ", 
      	              "argument 'packages' do not supply reasonable parametric", 
                      " distributions or some fatal error occured.\n",
              "To troubleshoot try:\n",
      	      "(1) changing argument input for 'packages'\n",
      	      "(2) adjusting rigorosity of 'sanity'\n",
      	      "(3) adjusting 'timeout'\n")
      	      df <- data.frame()
      	    } else {
              object <- sort(object, ic = ic)
              df <- data.frame(family = sapply(object@fits,
                                             function(f) f@family),
                             package = sapply(object@fits, function(f) 
                                                                  f@package),
                             ic = sapply(object@fits, function(f) f %@% ic),
                             params = sapply(object@fits, function(f) 
                               paste(names(f@estimatedValues),
                                     signif(f@estimatedValues,
                                            digits = 3),
                                     sep = " = ", collapse = '; ')))
              colnames(df) <- c("family", "package", ic, "params")
              df <- df[1:max(1,min(n, nrow(df))),]
      	    }
            return(new("globalfitSummary",
                       call = object@call,
                       data = object@data,
                       continuity = object@continuity,
                       method = object@method,
                       fits = df,
                       ic = ic
            ))
          }
)

setMethod(f = "show", signature = c("globalfitSummary"),
          def = function(object) {
            cat('\nCall: \n', object@call, '\n\n')
            if(is.null(object@continuity)) {
              cont <- ''
            } else if(object@continuity) {
              cont <-
                '\nAssumption: Data was generated from a continuous distribution.'
            } else if(!object@continuity) {
              cont <- 
                '\nAssumption: Data was generated from a discrete distribution.'
            }
            cat(length(object@data), 
                'data points entered. Distributions were fitted via', 
                object@method, 'estimation.',
                cont, '
                \nBest fits sorted by', object@ic, ':\n\n')
            print(object@fits, right = FALSE)
          }
)

setMethod(f = "print", signature = c("globalfitSummary"),
          def = function(x) {
            show(x)
          }
)

setMethod(f = "show", signature = c("globalfit"),
	  def = function(object) {
	    show(summary(object))
	  }
)


IC <- function(object, ic = "AIC", n = NULL) {
  if(is.null(n))
    n <- Inf
  if(!is.natural(n))
    stop("Argument 'count'  must be positive integer.")
  object <- sort(object, ic = ic)
  n <- min(length(object@fits), n)
  object@fits <- object@fits[1:n]
  x <- sapply(object@fits, function(object) object %@% ic)
  names(x) <- paste(sapply(object@fits, function(object) object@package), 
                    sapply(object@fits, function(object) object@family),
                    sep = "::")
  return(x)            
}

setMethod(f = "AIC", signature = c("globalfit"),
          def = function(object, n = Inf) {
            IC(object, ic = "AIC", n = n)
          })

setMethod(f = "BIC", signature = c("globalfit"),
          def = function(object, n = Inf) {
            IC(object, ic = "BIC", n = n)
          })

# setMethod(f = "AICc", signature = c("globalfit"),
#	  def = function(object, n = Inf) {
#	    IC(object, ic = "AICc", n = n)
#	  })

setMethod(f = "hist", signature = c("globalfit"),
          def = function(x, which = 1,  ic = c('BIC', 'AIC')) {
            ic <- match.arg(ic)
            if(is.null(which) || !is.numeric(which) ||
               abs(as.integer(which)) != which)
              stop("Argument 'which' must be a positive integer.")
            
            x <- sort(x, ic = ic)
            if (which > length(x@fits)) stop(
              "value of 'which' larger than the number of available results")
             
            lower <- min(x@data) - 0.2 * (max(x@data) - min(x@data))
            upper <- max(x@data) + 0.2 * (max(x@data) - min(x@data))

            selected_fit <- x@fits[[which]]
            if(x@continuity) {
              supporting_point <- seq(lower, upper, length.out = 300)
            } else {
              supporting_point <- seq(floor(lower), ceiling(upper))
            }
            breaks <- if (x@continuity) sqrt(length(x@data)) else 
                      0.5 + (min(x@data) - 1) : max(x@data)
            
            fun <- get_fun_from_package(type = "d", family = selected_fit)
            param_list <- split(selected_fit@estimatedValues, 
                                names(selected_fit@estimatedValues))
            param_list$x <- supporting_point
            density <- do.call(fun, param_list)
            
            # first create the histogram without plotting, because we need its 
            # maximum density for setting ylim
            h <- hist(x = x@data, breaks = breaks, plot = FALSE)
            
            plot(h, xlim = range(lower, upper), xlab = 'x', freq = FALSE,
                 ylim = range(0, max(h$density), max(density)), 
                 ylab = 'density',
                 main = paste0('Histogramm with density of \n',
                             selected_fit@package, '::', selected_fit@family)
            )
            
            # add true density (in discrete case add points 
            # at the possible x values)
            lines(supporting_point, density, col = 'green', lwd = 2)
            if (!x@continuity) points(supporting_point, density, 
                                      col = 'green', lwd = 2)
            
            h$estimation_points <- supporting_point
            h$estimated_density <- density
            invisible(h)
          }
        )
