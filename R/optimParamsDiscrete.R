## Authors 
## Moritz Kern, mkern@mail.uni-mannheim.de
## Benedikt Geier, bgeier@mail.uni-mannheim.de
## Borui N. Zhu, bzhu@mail.uni-mannheim.de
##
## Optimise the discrete and continuous parameters of a distribution
##
## Copyright (C) 2019 -- 2020 Moritz Kern, Benedikt Geier, Borui N. Zhu
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



### STATEMENT OF PURPOSE ###
############################

# This function is a wrapper for optimParamsContinuous
# Implementation is not final, target for now is to get a working solution 
# - no matter the runtime efficiency.
# This wrapper discerns the following cases:
# (a) no integer parameters: call optimParamsContinuous straight away
# (b) one integer parameter: use greedy solution. Questions to: Moritz Kern
# (c) multiple integer parameters: use naive solution. Question to: Borui N. Zhu

# TODO: in general, we often restrict parameters to be in a reasonable bound: 
#       [-100,100] (refer to case (b) )
# this might not hold. Maybe we can let the user specify something 
# like parscale, that lets us space out everything,
# and some kind of accuracy so we optimise rigorously locally 
# after a rough "global" view

# Function for doing the optimisation of discrete parameters
optimParamsDiscrete <- function(data, family, family_info, method = 'MLE',
                                prior = NULL, log = TRUE,
                                optim_method = 'L-BFGS-B', n_starting_points = 1,
                                debug_error = FALSE, show_optim_progress = FALSE,
                                on_error_use_best_result = TRUE, 
                                max_discrete_steps = 100, discrete_fast = TRUE,
                        				plot = FALSE, max_zoom_level = 4,
                        				timeout = 15) {
  
  # update defaults with priors
  if(length(prior) > 0) {
    match <- match(names(prior), names(family_info$lower))
    family_info$defaults[match] <- prior
  }

  #: define loop variables
  i <- 1
  
  # CASE 1: No discrete params -> we can directly redirect to 
  #optimParamsContinuous
  if (all(family_info$accepts_float)) {
    optim_res <- tryCatch({
      optimParamsContinuous(data = data, family = family, lower = family_info$lower,
                            upper = family_info$upper,
                            defaults = family_info$defaults, method = method,
                            fixed = c(), log = log, optim_method = optim_method,
                            n_starting_points = n_starting_points,
                            debug_error = debug_error,
                            show_optim_progress = show_optim_progress,
                            on_error_use_best_result= on_error_use_best_result,
                            timeout = timeout) 
      }, error = function(e) {
        message(e)
        return(NULL)
      }
    )
    if (is.null(optim_res)) return(NULL)
    
  } # CASE 2: exactly one discrete parameter 
  else if( sum(!family_info$accepts_float) == 1 ) {

    # get the discrete parameter
    dispar_id <- !family_info$accepts_float
    dispar_default <- family_info$defaults[dispar_id]
    
    # vector of the single optimisation results
    cont_optim_results <- vector('list',length = max_discrete_steps)       

    # history dataframe where the optim progress is stored
    history <- data.frame(matrix(NA, nrow = max_discrete_steps, ncol = 3))  
    colnames(history) <- c("param_value", "direction", "log_lik")
    
    # as we iterate both left and rightwards staring from the default value we 
    # need to save the current values
    # first iteration with the default value is considered to be left
    cur_left_val <- dispar_default
    cur_right_val <- dispar_default + 1
    # whether our left or right iteration has reached the border
    touched_lower <- touched_upper <- FALSE

    repeat {
      if (!touched_lower && ! touched_upper) {
        # make sure that first one is left!, 
        # apart from that always alternate if possible
        direction <- c("right", "left")[i%%2 + 1]  
      } else if (!touched_upper) {
        direction <- "right"
      } else if (!touched_lower) {
        direction <- "left"
      } else break
      
      dispar <- if (direction == "right") cur_right_val else cur_left_val
      if(show_optim_progress) cat("Discrete Parameter:", dispar, "\n")
      
      # optimise with the current fixed value
      curr_res <- tryCatch(optimParamsContinuous(
          data = data, family = family, lower = family_info$lower[!dispar_id],
          upper = family_info$upper[!dispar_id],
          defaults = family_info$defaults[!dispar_id], method = method,
          fixed = dispar, log = log, optim_method = optim_method,
          n_starting_points = n_starting_points, debug_error = debug_error,
          show_optim_progress = show_optim_progress,
          on_error_use_best_result = on_error_use_best_result,
	  timeout = timeout),
        error = function(e) {
          if(debug_error) message(e);
          NULL}
	)
      
      # if successful add results to dataframe
      if(!is.null(curr_res)) {
        cont_optim_results[[i]] <- curr_res
        cont_optim_results[[i]]$par[names(dispar_default)] <- dispar[1]
        history[i, ] <- list(dispar[1], 
                              direction, cont_optim_results[[i]]$value) 
      } else {
        history[i, ] <- list(dispar[1], direction, NA)
      }
      
      # update current iteration values and check whether bound is reached 
      # or score has not improved
      if (direction == "right") {
        cur_right_val <- cur_right_val + 1
        
        # get all the results for "right" achieved up to now
        # we stop when the current result is worse than 
        # the best one achieved in this direction
        relevant_hist <- history[history$direction == "right", "log_lik"]
        relevant_hist <- relevant_hist[!is.na(relevant_hist)]
        
        # stop when greater than upper limit or worse than best result in the
        # same direction
        touched_upper <- cur_right_val > family_info$upper[dispar_id] || (
          discrete_fast && length(relevant_hist) > 2 &&
            relevant_hist[length(relevant_hist)] < max(relevant_hist)
        )
      }
        
      # same for left direction  
      if (direction == "left") {
        cur_left_val <- cur_left_val - 1
        
        relevant_hist <- history[history$direction == "left", "log_lik"]
        relevant_hist <- relevant_hist[!is.na(relevant_hist)]
        
        touched_lower <- cur_left_val < family_info$lower[dispar_id] || (
          discrete_fast && length(relevant_hist) > 2 && 
            relevant_hist[length(relevant_hist)] < max(relevant_hist)
        )
      }
      
      i <- i+1
      
      # stop if not converged so far
      if(i > max_discrete_steps) {
        warning('Discrete Optimization aborted, did not converge.')
        break 
      }
    }
    
    if(show_optim_progress) print(history)
    if(plot) {
      plot(history$param_value, ifelse(is.finite(history$log_lik),
                                        history$log_lik, NA), 
           ylab = "log_lik", xlab = names(family_info$lower)[dispar_id])
    }
    
    # take the best result
    if (sum(!is.na(history$log_lik)) > 0) {
      optim_res <- cont_optim_results[[which.max(history$log_lik)]]
    } else {
      message("No valid discrete optimization result achieved")
      return(NULL)
    }
  } else { # Cases 3 & 4: more than one non-integer parameter
    non_floats <- !family_info$accepts_float
    num_discrete <- sum(non_floats)
    final_ll <- numeric(1)
    ## naive implementation
    # get parameter ranges of non-float parameters, make compact grid
    # To deal with a piori arbitrarily large values, 
    # let's try something I'll call the Google Earth algorithm
    # start with reasonable ranges
    zoom <- zoom_level <- rep(0, times = num_discrete)
    # at the start, centre over defaults. 
    # later: centre over maximum and zoom in/out
    centre <- family_info$defaults[non_floats]
    while_counter <- 0
    repeat {
      while_counter <- while_counter + 1
      zoom_level <- zoom_level + zoom
      # centre is always an integer
      ## Konstanten nicht mitten im Code
      grid_low <- centre-(25*(10^zoom_level)) 
      grid_high <- centre+(25*(10^zoom_level))
      stepsize <- rep(1, times = num_discrete)*(10^zoom_level)
      if(show_optim_progress) {
	cat('current zoom level:', zoom_level, '\n')
        cat('current focal point:\n')
        print(centre)
        cat('stepsizes:', stepsize, '\n')
      }
      lows <- pmax(family_info$lower[non_floats], grid_low)
      # at the start, centre over defaults. 
      # later: centre over maximum and zoom in/out
      highs <- pmin(family_info$upper[non_floats], grid_high)
      # get_params shall insure that lower and upper are all integers
      # is there a vectorised version of seq()?
      ## was macht hier Vectorize ?
      seq_vec <- Vectorize(seq.default, 
                           vectorize.args = c("from", "to", "by"), 
                           SIMPLIFY = FALSE)
      grid <- seq_vec(from = lows, to = highs, by = stepsize)
      grid <- expand.grid(grid)
      # output is a list, 
      # list entry number = position of param in family_info$lower 
      colnames(grid) <- names(family_info$lower)[non_floats]
      # number of columns of result matrix: 
      # number of variable parameters + two (loglikelihood & convergence code)
      num_free_params <- length(family_info$lower) - sum(non_floats)
      if(num_discrete < length(family_info$lower)) { # case 3
        grid_results <- matrix(NA, nrow = nrow(grid), ncol = 2 )
        colnames(grid_results) <- c("loglik", "convergence")
        ## pb <- txtProgressBar(min = 0, max = nrow(grid))
        ## print("Please stand by shortly...")
        for(i in 1:nrow(grid)) {
	  optim_res <- tryCatch(
	  {
	    optimParamsContinuous(data = data, family = family,
                                  lower = family_info$lower[!non_floats],
                                  upper = family_info$upper[!non_floats], 
                                  defaults = family_info$defaults[!non_floats],
                                  method = method, fixed = grid[i, ],
                                  prior = prior, log = log, 
                                  optim_method = optim_method,
                                  n_starting_points = n_starting_points,
                                  debug_error = debug_error, 
                                  show_optim_progress = show_optim_progress,
                                  on_error_use_best_result =
                                  on_error_use_best_result,
                                  no_second = TRUE, timeout = timeout)
	  },
	  error = function(e) {
	    message(e);
	    # generate a NA row of appropriate length to impute into grid_results
	    list(
	      val = NA,
	      convergence = 99
	    )
	  } # end error handler
	  ) # end tryCatch
	  grid_results[i, ] <- c(optim_res$val, optim_res$convergence)
	 ## setTxtProgressBar(pb, i)
  } # end for-loop in grid
        # drop all weird cases
        discrete_results <- grid_results[grid_results[,'convergence'] == 0, ]
      } else { # case 4
	grid_results <- matrix(NA, nrow = nrow(grid), ncol = 1)
        colnames(grid_results) <- c("loglik")
	##pb <- txtProgressBar(min = 0, max = nrow(grid))
	# print("Please stand by shortly...")
	for(i in 1:nrow(grid)) {
	  loglik_fun <- loglik(family = family, data = data, fixed = grid[i, ],
                               log = log, lower = family_info$lower,
                               upper = family_info$upper)
	  grid_results[i,] <- tryCatch(
	    {
	      loglik_fun()
	    },
	    error = function(e) {
	      NA
	    }
	  )
	 ## setTxtProgressBar(pb, i)
	}
      }
      optimum_index <- which.max(grid_results[,'loglik'])
      # cat('\noptimal index:', optimum_index, '\n') # for debug only
      final_ll <- grid_results[optimum_index,'loglik']

      # print(zoom_level)
      # print(max(zoom_level))
      
      # Google Earth: check if optimum is at the bound of our grid. 
      # If so, zoom out and center! if not: accept and break.
      optimum_gridcell <- grid[optimum_index, ]
      boundary_check <- ( (optimum_gridcell == grid_low) |
                          (optimum_gridcell == grid_high) )
      if (sum(boundary_check) > 0) {
	centre <- optimum_gridcell
        # zoom out only in dimensions where max was at boundaries
        zoom <- as.numeric(boundary_check)
	if(max(zoom_level) >= (max_zoom_level - 1) ) {
	  # par_outofbound <- non_floats[which.max(zoom_level)] 
	  # here which.max is okay, showing one parameter only is fine
	  # cat('Maximum zoom reached. Parameter', names(par_outofbound)[1], 
	  # 'appears to be far off.\nMaybe adjust priors or max_zoom_level?\n')
	  # Above is unstable as best parameters can be jointly off, and display 
	  # the name of a parameter that incidentally might be an okay one.
	
	  warning('Maximum zoom reached; some parameters appear to be far off.\n
	          Maybe adjust priors or max_zoom_level?\n')
	  break # no need to zoom in if par -> Inf
	} else {
	  if(show_optim_progress) print("Zooming out!")
	}
      } else if (max(zoom_level) > 0) {
	centre <- optimum_gridcell
        # since no boundary optima, zoom in wherever zoom is highest
	# do not use which.max, as index may not be unique
	which_max <- (zoom_level == max(zoom_level)) 
  zoom <- as.numeric(which_max) * (-1)
	if(show_optim_progress) print("Zooming back in!")
      } else { ## hmm. Unklar # BNZ: Der einzige Fall, der bis hierhin 
                              # überlebt, ist wenn die optimale Gitterzelle 
                              # nicht am Rande steht
	                      # UND der Zoom wieder auf genauester Ebene ist. Deswegen 
                        # ist hier das Optimum erreicht --> break
	break
      }
    } # while end

    # run optimParamsContinuous again for the best grid cell to retrieve 
    # information criteria, otherwise grid_results would blow up too much
    # difference to the optimParamsContinuous above: argument fixed is changed!
    if(num_discrete < length(family_info$lower)) {
      optim_res <- tryCatch(
        {
          optimParamsContinuous(data = data, family = family, 
                                lower = family_info$lower[!non_floats], 
                                upper = family_info$upper[!non_floats], 
                                defaults = family_info$defaults[!non_floats], 
                                method = method, 
                                fixed = grid[optimum_index, ], 
                                prior = prior, 
                                log = log, 
                                optim_method = optim_method, 
                                n_starting_points = n_starting_points, 
                                debug_error = debug_error, 
                                show_optim_progress = show_optim_progress, 
                                on_error_use_best_result =
			        on_error_use_best_result,
			        timeout = timeout)
        },
        # error should not occur because the combination 
        # had passed the first time!
        error = function(e) {
          message(e);
          list(par = rep(NA, times = (num_free_params)),
  	     val = NA,
  	     convergence = 99)
        }
      )
      final_params <- as.vector(c(optim_res$par, grid[optimum_index, ]), 
                                mode = "numeric")
      names(final_params) <- c(names(optim_res$par), 
                               colnames(grid[optimum_index, ]))
      reorder <- match(names(final_params), names(family_info$lower))
      optim_res$par <- final_params[reorder]
    } else {
      optim_res <- list()
      optim_res$par <- grid[optimum_index,]
    }
    optim_res$value <- final_ll
  }

  # ICs are the same, since discrete parameters are 
  # still parameters we optimise over
  ic <- informationCriteria(ll = optim_res$value, n = length(data),
                            k = length(family_info$upper))
  optim_res <- c(optim_res, ic)
  return(optim_res)
}
