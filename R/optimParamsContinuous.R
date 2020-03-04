## Authors 
## Moritz Kern, mkern@mail.uni-mannheim.de
## Benedikt Geier, bgeier@mail.uni-mannheim.de
##
## Optimise the continuous parameters of a distribution
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



# helper function, extracts the best result from all optimisation tries
get_best_result_from_progress <- function(optim_progress, param_names) {
  best_idx <- which.max(optim_progress$log_lik)
  best_row <- optim_progress[best_idx,]
  optim_result <- list()
  optim_result$value <- best_row$log_lik
  optim_result$par <- unlist(best_row[param_names])
  optim_result$convergence <- 51  # 51 indicates the warning code from optim
  
  return(optim_result)
}


## Parameters of optimParamsContinous:
# family: list with two elements "family" and "package"
# lower, upper and start_parameters must only contain 
# the continuous parameters that should be optimized
# prior: user-given prior information on parameters, 
# updates default values from get_param
# debug_error: show optimization progress when an error occured
# show_optim_progress: always show optimization progress
# on_error_use_best_result: if TRUE and an error occured during 
# optimization the best result achieved prior to the error will be taken
# n_starting_points: how many different starting points should 
# be used for optimisation. The best result will be taken.
optimParamsContinuous <- function(data, family, lower, upper, defaults,
                                  method = 'MLE', fixed = list(), prior = NULL,
                                  log = TRUE,
                                  optim_method = 'L-BFGS-B',
                                  n_starting_points = 1,
                                  debug_error = FALSE,
				  show_optim_progress = FALSE,
                                  on_error_use_best_result = TRUE,
                                  no_second = TRUE, timeout = 15) {

  if(method != 'MLE')
    stop('Not implemented.')
  if(length(lower) != length(upper) || length(defaults) != length(upper))
    stop('Length of lower and upper bounds vector do not coincide.')
  if(length(lower) == 0) {
    stop('No parameters to optimize as no bounds delivered.')
  }
  if(any(names(lower) != names(upper)) || 
     any(names(lower) != names(defaults)) ) {
    stop('Parameter names of lower and upper bounds and 
         start parameters must coincide. ')
  }
  
  stopifnot(n_starting_points >= 1)
  
  # replace default values from get_params with user-given priors
  if(length(prior) > 0) {
    prior_positions <- match(names(prior), names(lower), nomatch = NULL)
    # nomatch should not occur due to the check above: 
    # "Parameter names given as prior unknown"
    defaults[prior_positions] <- prior
  }
  
  # create dataframe where to save the optimization progress
  # 1 column for each parameter and a column for the associated log likelihood
  optim_progress <- data.frame(matrix(nrow = 0, 
                                      ncol = length(lower) + length(fixed) + 1))
  colnames(optim_progress) <- c(names(lower), names(fixed), "log_lik")
  
  on.exit({
    if ((show_optim_progress || (debug_error && !optim_successful))) {
      cat("Optimization progress:\n")
      print(tail(optim_progress, 2))
    }
  })
  
  # try multiple starting points hoping for a better result 
  # in case of multiple local minima
  optim_results <- vector("list", n_starting_points)

  for (i in 1:n_starting_points) {
    
      start_params <- if(i>1) sample_params(family, list(lower = lower, 
                                  upper = upper, accepts_float =! is.na(lower)), 
                                  params = lower) else defaults
      # cat("Sampling start parameters, Iteration:", i, "\n")
      # print(start_params)
    
      # Optimize first time
      if(show_optim_progress)
        cat("First Optimisation\n")
      
      # construct loglikelihood function, that only depends on the parameters
      loglik_fun <- loglik(family = family, data = data, fixed = fixed, 
                           log = log, upper = upper, lower = lower)
      safety_bound <- 1e-10
      
      if(is.numeric(timeout)) setTimeLimit(cpu = timeout, elapsed = timeout, 
                                           transient = TRUE)

      optim_result <- try(
        optim(start_params, loglik_fun, 
              control = list(fnscale = -1, trace = 0),
              lower = lower + safety_bound, 
              upper = upper - safety_bound, method = optim_method)
        , silent = TRUE)

      setTimeLimit(cpu = Inf, elapsed = Inf)
      
      if (is(optim_result, "try-error")) {
        optim_successful <- FALSE
        if(!on_error_use_best_result || nrow(optim_progress) == 0) {
          stop(optim_result, " occured during (first) optimisation for family ", 
               family$family,", fatal.\n")	
        } else {
          if (debug_error) message(optim_result, 
                  " occured during (first) optimization for family ",
                  family$family, 
                  " trying to take best result achieved up to now\n")
          # getting best result from optimization progress up to now
          optim_result <- get_best_result_from_progress(optim_progress, 
                                                param_names = names(lower))
        }
      } else {
        optim_successful <- TRUE
      }
      
      if(optim_result$convergence != 0 && debug_error) {
        print(tail(optim_progress, 2))
        warning('No convergence in first optimization for family ', 
                family$family)
      }
      
      
      ## NOTE: Second optimization has been removed for now
      # as parscale and fnscale have no positive effects on optimisation speed
      
      
      optim_results[[i]] <- optim_result
  } # end for-loop over starting values
  
  # extract best optim result of all the starting points
  best_idx <- which.max(sapply(optim_results, function(x) x$value))
  optim_result <- optim_results[[best_idx]]
  
  if (nrow(optim_progress) > 0 &&
      optim_result$value < max(optim_progress$log_lik, na.rm = TRUE) - 1e-8) {
    if(debug_error) {
      message("Final Optimization result for family ", family$family, 
            " is worse than the best result achieved during optimization")
      cat("Diff to best:", abs(optim_result$value - max(optim_progress$log_lik)),
          "\n")
    }
  }
  
  return(list(
    par = optim_result$par,
    value = optim_result$value,
    convergence = optim_result$convergence
    )
  )
}


# TODO: set fnscale and parscale appropriately -> DEACTIVATED right now

# TODO: optional and ADDITIONAL transformation of data for families with
#       bounded support, such as stats::beta or stats::unif

# NOTE: optim_progress contains all calls to log_lik-function, 
#       that is also the calls for estimating the gradient, thats why usually
#       always 5 rows refer to the same optimization step

