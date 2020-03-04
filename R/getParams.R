## Authors 
## Benedikt Geier, bgeier@mail.uni-mannheim.de
##
## Get properties of a single distribution family and its parameters
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



# ----------------------------------------------------------------------
# 1) Try to get the parameters (+ infos) from a distribution
# ----------------------------------------------------------------------

## Main ideas:
# 1) extract all possible params from r... function as the r function should 
#    need all params
# 2) check if all of them have default values set
# 2.1) If yes go to 3)
# 2.2) If not for each param with missing default value guess such a value 
#     and test if r..(1, params) returns a value or NA
#      If a valid value is returned take the current set of default values, 
#      otherwise try a different combination of default values
# 3) For each of the params test some values (non-integer, negative, 
#    not in [0,1]) while keeping the others at their defaults

## optional TODO:
# check whether param name contains prob or sth like that -> range should be[0,1]

# fam <- list(package = "stats", family = "beta") # gamma

# ----------------------------------------------------------------------
# (1.1) Given distribution family, return list of parameters
# ----------------------------------------------------------------------

get_all_params <- function(fam) {
  # idea: all params need to be present in the r... method for generating random
  # samples from the distribution
  
  fun <- get_fun_from_package(type = "r", family = fam)
  all_params <- formals(fun)
  
  # if "nn" is contained then "n" is probably a real param , c.f. rhyper
  to_remove <- if("nn" %in% names(all_params)) c("", "nn") else c("", "n")
  
  # remove empty and the n_samples argument
  all_params <- all_params[! names(all_params) %in% to_remove]
  
  # in cases like gamma distribution we have multiple parameters 
  # describing the same (e.g. rate and scale)
  # -> we need to drop those as we can't optimize them independently
  for (i in 1:length(names(all_params))) {
    param <- names(all_params[i])
    # TRUE if default val is function of another param, e.g. scale = 1/rate
    if (typeof(all_params[[param]]) == "language") {        
     # divide the function into its components
     function_parts <- as.character(all_params[[param]])      
     
     # check if one function part is the name of one of the earlier parameters,
     # and remove if thats true
     if (i > 1 && any(function_parts %in% names(all_params)[1:(i-1)]))
       cat(param, "is realated to a parameter listed before and 
           thus will be removed.\n")
       all_params[[param]] <- NULL
    }
  }
  return(all_params)
}

# data structure "all_params"
# list, for given distribution family;
# list elements name: name of parameter
# list elements field: default value, NULL if not set

# -----------------------------------------------------------------------------
# (1.2) Given distribution family, parameters, 
#       some x-values: test if combination is valid
# -----------------------------------------------------------------------------

validate_values <- function(fam, n_or_nn, params, x_test) {
  # try to generate a random number from the distribution
  rfun <- get_fun_from_package(type = "r", family = fam)
  dfun <- get_fun_from_package(type = "d", family = fam)
  
  r <- do.call(rfun, c(n_or_nn, params))
  # additionally check whether values are valid for density function
  # as this sometimes takes a while we use a timeout to stop the execution 
  # after a while
  # however we consider a timeout as a valid parameter value as no error 
  # occurs (it just takes too long)
  if (is.finite(r)) {
    r_ <- eval_with_timeout(do.call(dfun, c(x_test, params)), timeout = 1, 
                            return_value_on_timeout = "TIMEOUT")
    
    if (any(r_ == "TIMEOUT")) message(fam, " produced timeout for params ", 
                                 paste(names(params), params, sep = ": ", 
                                 collapse = ","), " on ", x_test)
    return(any(!is.na(r_)))
  } else {
    return(FALSE)
  }
}

# return value: boolean
# does parameter combination in input yield a valid value?

# -----------------------------------------------------------------------------
# (1.3) Given updated list of parameters, and a specific family, 
#       find one parameter combination that works
# -----------------------------------------------------------------------------

# first we need to find one set of values for each of the params that actually 
# works before we can check for valid values for each
# of the single params individually

get_default_values <- function(all_params, fam) {
  # missing values seem to have type "symbol"
  missing_defaults <- sapply(all_params, function(x) typeof(x) == "symbol")   
  
  if (sum(missing_defaults) == 0) return(all_params)
  
  with_defaults <- all_params[!missing_defaults]
  non_defaults <- all_params[missing_defaults]
  
  # both floats and integers and both positive and negative as well as 0 so 
  # that at least one of those hopefully is valid
  # start with 0.5 to avoid extreme values in case of probabilities
  default_guesses <- c(0.5, 1, 0, -0.5, -1)
  
  # create a dataframe with all combinations of default guesses
  combs <- expand.grid(lapply(non_defaults, function(x) default_guesses))
  combs_list <- split(combs, seq(nrow(combs)))   # convert to list for iteration
  valid_params <- NULL
  
  # parameter that describes the number of random numbers to take, 
  # usually "n", but in cases like hyper "nn"
  rfun <- get_fun_from_package(type = "r", family = fam)
  n_or_nn <- if (! "nn" %in% names(formals(rfun))) list(n = 1) else list(nn = 1)
  x_test <- list(x = seq(-10, 10, 1))
  
  errors <- c()
  for (i in 1:length(combs_list)) {
    
    # combine fixed with guessed default values 
    # and try to generate a random number
    # in most cases invalid parameter choices will just generate a warning 
    # and return NA, but sometimes also an error is thrown
    # so we need to handle both
    curr_params <- c(with_defaults, combs_list[[i]])
    res <- suppressWarnings(tryCatch({
      validate_values(fam, n_or_nn, curr_params, x_test)
    },
    error = function(e) {
      errors <<- union(errors, strsplit(as.character(e), ":", 
                                        fixed = TRUE)[[1]][2])
      return(FALSE)
    }))
    
    # break if we've found a set of valid values
    if (res){
      valid_params <- curr_params
      #cat("Found the following set of valid default values 
      # for family", fam$family, ":", 
      #    paste(names(valid_params), valid_params , sep=": ", 
      #          collapse=", "), "\n")
      break
    } 
  }
  if (is.null(valid_params)) {
    
    message("Could not find a set of valid default values for family ", fam,
            "\nErrors:", errors)
  }
  return(valid_params)
}
# return value: list with names corresponding to the parameter names and 
# values to their default values


#------------------------------------------------------------------------------ 
# (1.4) Get parameter ranges
# -----------------------------------------------------------------------------

# Now we can iterate over the parameters while keeping fixed all others 
# in order to guess some valid ranges
# Note that all_params will now always contain some default values for each of 
# the params (as long as one has been found)

# ----------------------------------------------------------------------------- 
# function that checks if all "values" are valid for "param". "all_params" 
# contains the values of the other params of the family "fam"
# ----------------------------------------------------------------------------- 

check_values_for_param <- function(param, all_params, fam, values) {
  
  # parameter that describes the number of random numbers to take, 
  # usually "n", but in cases like hyper "nn"
  rfun <- get_fun_from_package(type = "r", family = fam)
  n_or_nn <- if (! "nn" %in% names(formals(rfun))) list(n = 1) else list(nn = 1)
  x_test <- list(x = seq(-10, 10, 1))
  
  res <- suppressWarnings(
    sapply(values, function(x) {
      # overwrite the value of "param" to each of the "values" 
      # that should be tested
      all_params[[param]] <- x;
      tryCatch(
        validate_values(fam, n_or_nn, all_params, x_test),
      error = function(e) return(FALSE))
    })
  )
  return(res)
}
# return value: TRUE if value was valid, FALSE if not

# --------------------------------------------------------------------------- 
# function that iterates over descending step sizes to find the minimal and 
# maximal valid value of a parameter
# --------------------------------------------------------------------------- 

# Explanation:
  # in previous iteration we tested the values [-10, -5, 0, 5, 10], 
  # that is step_size was 5
  # valid values: [FALSE, FALSE, TRUE, TRUE, TRUE]
  # now the lower limit can be anywhere between -5 and 0
  # so we test with next step size=1 the values [-5, -4, -3, -2, -1, 0]
  # from those values we again search the minimum valid value and continue
  # in the same way with always smaller steps
  # for the upper bound the procedure is the same, just that we search for 
  # valid values in the range higher than the current upper limit

iterate_min_max_vals <- function(param, all_params, fam, cur_val, 
                                 step_sizes, is_min = TRUE) {
  for (i in 1:(length(step_sizes)-1)) {
    
    # first get the interval that should be tested (lower than cur_val when 
    # searching lower limit and higher than cur_val when searchin upper limit)
    if (is_min) {
      border <- cur_val - step_sizes[i]
      vals <- seq(border, cur_val, by = step_sizes[i+1])
    } else{
      border <- cur_val + step_sizes[i]
      vals <- seq(cur_val, border, by = step_sizes[i+1])
    }
    # test values and adjust the current estimate
    # cat(param, "-> Currently checking values in the interval 
    # [", min(vals), ",", max(vals), "] with step size", step_sizes[i+1], "\n")
    check_res <- check_values_for_param(param, all_params, fam, vals)
    
    # usually at least the first (when is_min = FALSE) 
    # or the last (when is_min = TRUE) should be valid
    # but due to overflow it can happen that none is TRUE 
    # then we break with the current value
    if (!any(check_res)) return(cur_val)
    
    cur_val <- if(is_min) min(vals[check_res]) else max(vals[check_res])
    
    # if the smallest (or highest) of the tested values was valid we can break
    # however this should usually not happen when the method is used as below 
    # (but happens due to float imprecisions...)
    if (cur_val == border) {
      # cat("Abbruchbedingung")
      return(cur_val)
    }
  }
  # cat(param, "-> Final value:", cur_val, "\n")
  return(cur_val)
}


# --------------------------------------------------------------------------- 
# (1.5) Main function for generating the info for each of the params
# ---------------------------------------------------------------------------

get_param_ranges <- function(all_params, fam) {
  
  # SPECIAL CASE: uniform distribution
  if (fam$family == "unif") {
    lower <- rep(-Inf, length(all_params))
    upper <- rep(Inf, length(all_params))
    accepts_float <- rep(TRUE, length(all_params))
    names(lower) <- names(upper) <- names(accepts_float) <- names(all_params)
    
    return(list(lower = lower,
                upper = upper,
                accepts_float = accepts_float, 
                defaults = unlist(all_params)))
  }
  
  # create empty result vectors
  lower <- upper <- accepts_float <- rep(NA, length(all_params))
  names(lower) <- names(upper) <- names(accepts_float) <- names(all_params)
  
  # Parameters for the iteration for finding the valid ranges of each parameter
  initial_min_val <- -1e3
  initial_max_val <- 1e3
  step_sizes <- c(1e2, 50, 10, 1, 0.1, 0.01)
  n <- 10
  
  for (param in names(all_params)){
    
    # 1) Try to find upper and lower bounds (if there are some)
    
    # add the default value to have at least one valid entry
    vals <- seq(initial_min_val, initial_max_val, by = step_sizes[1]) + 
                all_params[[param]]  
    # cat("current step size:", step_sizes[1], "\n")
    check_res <- check_values_for_param(param, all_params, fam, vals)
    
    # if lowest or highest value was valid in the first check we already have a
    # min_val or max_val
    # otherwise we iterate with the above method
    min_val <- if (check_res[1]) -Inf else 
      {iterate_min_max_vals(param = param, all_params = all_params, fam = fam,
                            cur_val = min(vals[check_res]), 
                            step_sizes = step_sizes, is_min = TRUE)}


    max_val <- if (check_res[length(check_res)]) Inf else 
      {iterate_min_max_vals(param = param, all_params = all_params, fam = fam,
                            cur_val = max(vals[check_res]),
                            step_sizes = step_sizes, is_min = FALSE)}
    # set the estimated values in the named vactor that will be returned
    lower[param] <- min_val
    upper[param] <- max_val
    
    # 2) Check if the parameter accepts floats or only integers
    testsequence <- runif(n, max(min_val, -1e3), min(max_val, 1e3))
    testsequence <- unique(c(testsequence, trunc(testsequence)))
    num_tests <- length(testsequence)

    is_integer <- testsequence %% 1 == 0
    num_integer <- sum(is_integer)

    testoutcome <- check_values_for_param(param, all_params, fam, testsequence)
    accepted_int <- sum(is_integer & testoutcome)
    accepted_int_rate <- accepted_int/num_integer
    accepted_float <- sum(!is_integer & testoutcome)
    accepted_float_rate <- accepted_float/(num_tests - num_integer)

    accepts_float[param] <- accepted_float_rate > 1/(num_tests - num_integer)
    if (!accepts_float[param] && accepted_int_rate == 0) {
      stop("distribution does not seem to accept any values")
    }
  }
  return(list(lower = lower,
              upper = upper,
              accepts_float = accepts_float, 
              defaults = unlist(all_params)))
}

# return value: 4 component list
# each list entry: vector with length = number of parameters
#                  each list entry saves either lower, upper, or etc.


### Function that checks if log is working ------------------------------------
check_log <- function(fam) {
  dfun <- get_fun_from_package(type = "d", family = fam)
  return('log' %in% names(formals(dfun)))
}


### Function that checks whether a family is a discrete distribution, 
### that is it only takes integers as values ----------
check_integer <- function(fam, all_params) {
  rfun <- get_fun_from_package(type = "r", family = fam)
  n_test <- 10
  n_or_nn <-
    if (! "nn" %in% names(formals(rfun))) list(n = n_test)
    else list(nn = n_test)
  args_ <- c(all_params, n_or_nn)
  res <- do.call(rfun, args = args_)
  return(all(abs(res %% 1) < sqrt(.Machine$double.eps)))
}

# -----------------------------------------------------------------------------
# (2) Determining the support of a distribution
# -----------------------------------------------------------------------------

# Function for determining the support of a given distribution family and 
# whether the limits of the support are determined by one
# of the distributions parameters
# params should be a list containing named vectors lower, upper, accepts_float 
# with one entry for each of the distributions parameters
# additionally params$discrete specifies whether fam is a discrete distribution
get_support <- function(fam, params) {
  
  # Parameters used below
  cap_at <- 10
  ignore_extreme_vals_perc <- 0.1
  max_test_range <- 100           # range of x values to check the density on
  precision <- 0.01               # how exact should the support be estimated
  n_test <- 11                    # number of values to test for each parameter
  support_to_inf_limit <- 50
  
  # in case the parameter is unbounded we cap it to avoid extreme values
  low_capped <- pmax(params$lower, -cap_at)
  upp_capped <- pmin(params$upper, cap_at)
  
  # we additionally ignore the 10% highest and lowest parameter values 
  # (e.g. for "binom" we only consider prop in [0.1, 0.9])
  low <- ifelse(params$lower == low_capped, 
                low_capped + ignore_extreme_vals_perc * (upp_capped-low_capped), 
                low_capped)
  upp <- ifelse(params$upper == upp_capped, 
                upp_capped - ignore_extreme_vals_perc * (upp_capped-low_capped), 
                upp_capped)
  
  # initialize named vector that stores whether each parameter 
  # determines the bounds of a distribution
  supp_min_depends_on <- supp_max_depends_on<- rep(FALSE, length(params$lower))
  names(supp_max_depends_on) <- 
    names(supp_min_depends_on) <- names(params$lower)
  
  # define the base choices for all parameters that are chosen when only 
  # varying one parameter and keeping the others constant
  base_choices <- as.list(ifelse(params$accepts_float, params$defaults,
                                 round(params$defaults)))
  
  # define the sequence of test points
  x <- seq(-max_test_range, max_test_range, precision)
  if (params$discrete) x <- unique(round(x))
  
  # and add the points to the list of parameter choices
  base_choices$x <- x
  
  # initialize limits of support to maximum / minimum possible value each
  support_min <- Inf
  support_max <- -Inf
  
  dfun <- get_fun_from_package(type = "d", family = fam)
  
  for (param in names(low)) {
    # define n_test equally distributed values for the current param 
    # dependend on its adapted range from above
    param_choices <- seq(low[param], upp[param], length.out = n_test)
    
    if (!params$accepts_float[param]) param_choices <- trunc(param_choices)
    
    # copy base choices to args_ so that we can change 
    # the value for the current param below
    args <- base_choices
    
    get_result_mat <- function(param_choices){
      
      # row i of the result matrix will be the density values at x 
      # when taking the i-th choice for the current param
      result_mat <- matrix(NA, nrow = length(param_choices), ncol = length(x))
      for (i in 1:length(param_choices)) { 
        choice <- param_choices[i]
        
        # calulate density value and add to result matrix
        args[[param]] <- choice
        res <- suppressWarnings(do.call(dfun, args = args))
        result_mat[i, ] <- res
        i <- i+1
      }
      return(result_mat)
    }
    
    result_mat <- get_result_mat(param_choices)
    
    # for each row calculate the minimum and maximum evaluation point 
    # with positive density
    row_support_min <- apply(result_mat, 1, 
      function(row) {if (length(which(row > 0)) > 0) 
                        x[min(which(row > 0))] else Inf})
    row_support_max <- apply(result_mat, 1, 
      function(row) {if (length(which(row > 0)) > 0) 
                          x[max(which(row > 0))] else -Inf})
    
    # check if the lower or upper bound is always the same as the current 
    # parameter value (up to the chosen precision + some small machine error)
    # then the support depends on the current parameter and the support is at 
    # least as big as the possible ranges of this parameter
    # notice that we need to ignore rows where the min or max support value 
    # was +/- Inf
    min_criterium <- abs( param_choices - row_support_min )
    max_criterium <- abs( param_choices - row_support_max )
    if(max(min_criterium[is.finite(row_support_min)]) <= precision + 1e-10) {
      supp_min_depends_on[param] <- TRUE
      support_min <- min(params$lower[param], support_min)
      support_max <- max(max(row_support_max), support_max)
    }
    if(max(max_criterium[is.finite(row_support_min)]) <= precision + 1e-10) {
      supp_max_depends_on[param] <- TRUE
      support_max <- max(params$upper[param], support_max)
      support_min <- min(min(row_support_min), support_min)
    }
    
    # else we just adapt the current maximum and minimum support values with
    # the minimum or maximum row support
    if (!supp_min_depends_on[param] && ! supp_max_depends_on[param]) {
      support_min <- min(min(row_support_min), support_min)
      support_max <- max(max(row_support_max), support_max)
    }
    
    # cat("After param", param, "--> \tsupport_min:", support_min, 
    # "\tsupport_max", support_max, "\n")
  }

  # if minimum / maximum support is small / high enough we assume that the
  # support is the whole real line
  if (support_min <= -support_to_inf_limit) support_min <- -Inf
  if (support_max >= support_to_inf_limit) support_max <- Inf
  
  return(list(support_min = support_min, support_max = support_max, 
              supp_max_depends_on = supp_max_depends_on,
              supp_min_depends_on = supp_min_depends_on))
}

# -----------------------------------------------------------------------------
# (3) Final function
# -----------------------------------------------------------------------------

standardizeFam <- function(fam, package){ 
  if (missing(package) || length(package) == 0) {
    if (!is(fam, "optimParams") && !is.list(fam)) {
      if (!is.character(fam))
        stop("'fam' must be a character string or a list.")
      if (length(fam) == 2) {
        fam <- list(fam["family"], package = fam["package"])
      }
      if (length(fam) == 1) {
        names <- sapply(FamilyList, function(x) x$family)
        idx <- pmatch(fam, names)
        if (is.na(idx))
          stop("The family can not be identified without the explicitly 
               given argument 'package'")
        fam <- list(family = fam, package = FamilyList[[idx]]$package)
      } else stop("The length of 'fam' must be one")
    }
  } else {
    fam <- list(family = fam, package = package)
  }
  return(fam)
}

getParams <- function(fam, package){
  
  fam <- standardizeFam(fam, package)

  # 1) Get list of all parameters:
  all_params <- get_all_params(fam)
  
  # 2) Add default values to all params that don't have any
  all_params <- get_default_values(all_params, fam)
  
  if (is.null(all_params)) return(NULL)
  
  # 3) Get valid parameter ranges:
  result <- get_param_ranges(all_params, fam)
  
  # 4) Add log argument
  result$log <- check_log(fam)
  
  # 5) Add discrete argument
  result$discrete <- check_integer(fam, all_params)
  
  # 6) Add support informations
  supp <- get_support(fam , result)
  result <- c(result, supp)
  
  return(result)
}


### Open problems & TODO:
# 1) Distributions like nbinom where 2 params ("prob" and "mu") describe 
#    the same but only one may be set and 
#    none of them has a default value derived from the other
# -> IGNORED at the moment -> will be filtered out
# 2) Distributions like "unif" where the parameters interact 
#    -> ranges can be represented as [lower, upper] but rather as min <= max
#  -> SOLVED
# 3) Distribution hyper: here rhyper also works with floats for m,n,k 
#    but dhyper not, maybe also check d function in check_values_for_param
#    Current errors have to be catched in get_support but it would be better 
#    if they didn't occur at all
#    -> SOLVED
# 4) when scaling GLOBALFIT needs to consider whether support_max or support_min
#   depends on certain parameters, then it needs to adapt the upper and
#   lower bounds of those parameters using the data, 
#   i.e. for unif set upper["min"] <- min(data) before giving to optim_param
#   -> NOT DONE at the moment
# 5) extend to and test with other packages
#    -> DONE
