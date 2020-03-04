## Authors 
## Niclas Lietzow, nlietzow@mail.uni-mannheim.de
## Till Freihaut, tfreihau@mail.uni-mannheim.de
## Leonardo Vela, lvela@mail.uni-mannheim.de
##
## Define the log-likelihood function for given parameters
##
## Copyright (C) 2019 -- 2020 Niclas Lietzow, Till Freihaut and Leonardo Vela
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


loglik <- function(family, data, fixed = list(), log, lower, upper) {
  
  stopifnot(length(upper)>0)
  
  arguments <- list(x=data) 
  
  # check wheter log-distribution function is directly available for 
  # distribution
  if(log)
    arguments$log <- TRUE
  # add fixed parameter values of distribution to list
  if(length(fixed)>0) {
    arguments <- c(arguments, fixed)
  }
  
  
  # define loglikelihood function
  # BNZ: allow for empty params for distributions with all-integer parameters
  likelihood <- function(params = NULL) { 

    for(param_name in names(params)) {
      if(lower[param_name] > params[[param_name]] ||
         upper[param_name] < params[[param_name]])
         stop('Parameter ', param_name, ' with value ', params[[param_name]],
              ' outside the boundaries.')
    }
    
    # Add params with names of parameters to arguments list
    arguments <- c(arguments, params)

    summands <- do.call(get_fun_from_package(type = "d", family = family),
                        args = arguments)
    
    if(any(is.na(summands))) stop('In Log-Likelihood-Function NA occured.')
    
    # log values if not log so far
    if(!log) {
      summands <- log(summands)
      ## warning('Could be numerically instable.')
    } 
    loglik_value <- sum(summands)
    
    ## The following is only for tracking the optimisation progress 
    ## (might be deactivated sometime)
    # recursively go through parent frames and check whether there is a variable 
    # tracks the optimisation process

    # if yes then add a new row to the progress dataframe in the closest parent 
    # frame
    for (i in 1:length(sys.parents())) {
       
      if (exists("optim_progress", envir = parent.frame(i))) {
        # cat("Found optimization progress in parent frame", i, "\n")
        progress <- get("optim_progress", envir = parent.frame(i))
        #print(progress)
        #print(c(param_values, fixed, log_lik = ll))
        progress[nrow(progress)+1, ] <- c(params, fixed, log_lik = loglik_value)
        assign("optim_progress", envir = parent.frame(i), progress)
        # print(tail(progress,2))
        break
      }
    }
    
    # message(loglik_value)
    return(loglik_value)
  }
  return(likelihood)
}
