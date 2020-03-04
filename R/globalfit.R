## Authors 
## Moritz Lauff, mlauff@mail.uni-mannheim.de
## Kiril Dik, kdik@mail.uni-mannheim.de
## Moritz Kern, mkern@mail.uni-mannheim.de
## Nadine Tampe, ntampe@mail.uni-mannheim.de
## Borui Niklas Zhu, bzhu@mail.uni-mannheim.de
## Benedikt Geier, bgeier@mail.uni-mannheim.de
## Helene Peter, hpeter@mail.uni-mannheim.de
##
## Fit multiple distribution families to a given univariate dataset
##
## Copyright (C) 2019 -- 2020 
## Moritz Lauff, Kiril Dik, Moritz Kern, Nadine Tampe, Borui Niklas Zhu, 
## Benedikt Geier
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




### 1)  necessary helper functions --------------------------------------------

# determine decimals
getDecimals <- function(x){
  
  return(round(x %% 1, 10))
  
}


some_percent <- function(decs, numbers, percent){
  # some_percent: TRUE, if for a given number of decimals
  #               at least a proportion of 'percent' out of
  #               all possible decimals occurs 
  
  for (i in min(numbers):max(numbers)){
    if (length(unique(decs[numbers == i])) >= percent * 10){
      return(TRUE) 
    }
  }
  return(FALSE)
}


# test, if data is discrete
is.discrete <- function(data, border = 0.35, percent = 0.8){
  
  # convert data, if not a vector
  if(is.data.frame(data)){
    data <- as.vector(data[,1])
  }
  
  # remove NA's
  data <- data[!is.na(data)]
  
  obs <- length(data)
  decs <- getDecimals(data)
  numbers <- nchar(as.character(decs)) - 2
  numbers[numbers == -1] <- 1
  n_unique_dec <- length(unique(decs))
  
  if (1 / border > obs){
    border <- 1 / obs
  }
  
  if (n_unique_dec / obs <= border && 
      !any(numbers >= 4) && 
      !some_percent(decs, numbers, percent)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# treatment of discrete non-integral numbers
disc_trafo <- function(data){
  
  # transformation only if data is discrete
  if (is.discrete(data)){ 
    
    data_new <- sort(data) # sort the data
    # data with corresponding decimals
    data_new <- list(data_new = data_new, decimals = getDecimals(data_new))
    
    # occuring distinct decimals
    unique_decimals <- sort(unique(data_new$decimals))
    m <- length(unique_decimals) # number of occuring distinct decimals
    
    # decimals and transformed decimals
    decimals <- list(original_decimals = unique_decimals,
                     new_decimals = round(seq(0, (m-1)/m, 1/m), 10))
    
    # put together
    data_new <- merge(data_new, decimals,  
                      by.x = "decimals", by.y = "original_decimals") 
    
    data_new <- (data_new$data_new - data_new$decimals 
                 + data_new$new_decimals) * m # new data
    
    return(list(data = data_new,
                trafo_df = decimals,
                discrete = TRUE,
                trafo_decription = 
                  paste0("Divide the simulated data by ",
                         m,
                         " and replace the decimals c(", 
                         paste(decimals$original_decimals, collapse = ", "),
                         ") of the simulated data by c(",
                         paste(decimals$new_decimals, collapse = ", "),
                         ").")))
    
  } else {
    return(list(data = data,
                trafo_df = NULL,
                discrete = FALSE,
                trafo_decription = NULL))
  }
  
}


### 2) Main Function ----------------------------------------------------------

globalfit <- function(data, continuity = NULL, method = "MLE", verbose = TRUE,
                      packages = "stats", append_packages = FALSE,
                      cores = NULL, max_dim_discrete = Inf,
                      sanity = 1, timeout = 5) {
  
  # set debug to TRUE to get a log file containing the messages emitted during 
  # optimization
  debug <- FALSE
  
  ic <- "BIC"

  all_funs <- c('%@%', 'check_integer', 'check_log', 'check_values_for_param', 
                'construct_package_list', 'disc_trafo', 'eval_with_timeout', 
                'fitting_sanity_check', 'get_all_params', 
                'get_best_result_from_progress', 'get_default_values', 
                'get_fun_from_package', 'get_fun_from_package_internal', 
                'get_param_ranges', 'get_support', 'getDecimals', 'getFamilies', 
                'getFamily', 'getParams', 'globalfit', 'IC', 
                'informationCriteria', 'is.discrete', 'is.natural', 
                'iterate_min_max_vals', 'iterate_packages', 'loglik', 
                'optimParamsContinuous', 'optimParamsDiscrete', 'print', 
                'sample_data', 'sample_params', 'some_percent', 'sort', 
                'standardizeFam', 'validate_values', 'write_file')
  
  ## Input Validation

  if ( length(sanity) != 1 || !( (is.numeric(sanity) && sanity >= 0) ||
			 isFALSE(sanity)) ) {
    stop("Invalid input for argument 'sanity'.")
  }
  do_sanity <- !isFALSE(sanity)
  
  if( length(timeout) != 1 || !(is.logical(timeout) || is.numeric(timeout) ) ||
     (is.numeric(timeout) && timeout < 0) ||
     (is.logical(timeout) && timeout == TRUE))
    stop("Invalid input for argument 'timeout'")

  families <- FamilyList
  installed <- rownames(installed.packages())

  if(length(packages) > 0) {
    
    if(is.vector(packages) && typeof(packages) == "character") {
      
      missing_pkgs <- setdiff(packages, installed)
      if (length(missing_pkgs) > 0) {
        message("The following packages were provided to argument 'packages' 
                but are not installed, so they will be ignored. ",
                "Please install manually: ",
                paste(missing_pkgs, collapse = ", "))
      }
      
      packages <- intersect(packages, installed)
      
      known_packages <- unique(sapply(families, function(x) x$package))
      
      additionals <- setdiff(packages, known_packages)
      if (length(additionals) > 0) {
        message("The following packages were provided in argument 'packages' 
                but are not part of the default set of packages: ",
                paste(additionals, collapse=", "),
                "\nThus the distribution families in those packages need to be 
                extracted now, which might take some time. ",
                "When executed multiple times, consider extracting those 
                families once with 'getFamilies(packages)' ",
                "and provide the result of that to argument 'packages'.")
        additionals_info <- iterate_packages(additionals)
        if (length(additionals_info) == 0) {
          message("No distribution families found in the additionally 
                  provided packages.")
        }
      } else {
        additionals_info <- list()
      }
      
      # add the manual ones to FamilyList as used in default
      if(append_packages) {
        families <- c(families, additionals_info)
        
        # ignore whatever else is in FamilyList
      } else {
        # these are specified by the user, but params are known
        known <- packages[! packages %in% additionals] 
        known <- families[ which(sapply(families, 
                                        function(x) x$package %in% known)) ]
        families <- c(known, additionals_info)
      }
      
    } else if(is.list(packages)) {
      # in the future this should be deprecated in favour of an S4 object
      # with better validity check 
      if(append_packages) {
        families <- c(families, additionals)
      } else {
        families <- packages
      }
    } else {
      stop("Invalid argument 'packages'.")
    }
  }
  
  if (length(families) == 0) {
    stop("The provided input to argument 'packages' did not 
         contain any distribution family. Can not optimize.")
  }
  
  # filter out those distributions that have too many discrete parameters.
  if(max_dim_discrete < Inf) {
    families <- families[which(sapply(families, 
                              function(x) {
                                sum(x$family_info$discrete) <= max_dim_discrete
                                } )) ]
  }
  
  # indices for the discrete distributions
  discrete_families <- which(sapply(families, 
                                    function(x) x$family_info$discrete))
  
  # if not specified, find out whether data is continuous or not
  if (is.null(continuity)){
    trafo_list <- disc_trafo(data)
    data <- trafo_list$data
    continuity <- !trafo_list$discrete
  } 
  
  # the following simpler version is not working as vector[-integer(0)] = empty
  # relevant_families <-  families[if (continuity) -discrete_families else discrete_families]
  if (continuity) {
    relevant_families <- 
      if (length(discrete_families) > 0) families[-discrete_families] else families
  } else {
    if (length(discrete_families) > 0)
      relevant_families <- families[discrete_families]
    else
      stop("Data is assumed to be discrete but there are no discrete distributions", 
           "contained in the set of distributions to be compared")
  }
  
  # Again check that all of the families that should 
  # be compared are also installed
  all_pkgs <- sapply(relevant_families, function(x) x$package)
  all_pkgs_unique <- unique(all_pkgs)
  missing_pkgs <- setdiff(all_pkgs_unique, installed)
  if (length(missing_pkgs) > 0) {
    message("The following packages are not installed, and are thus ignored",  
            "during optimisation. If you want to use them please install",
	    "manually: ", paste(missing_pkgs, collapse=", "))
    relevant_families <- relevant_families[!(all_pkgs %in% missing_pkgs)]
  }
  
  if(length(relevant_families) < 1)
    stop("No relevant distribution families to fit or compare.")

  if(verbose)
    message("Comparing the following distribution families: ", 
            paste(sapply(relevant_families, function(x) x$family), 
                  collapse = ", "))
  
  if(is.null(cores)) {
    CRAN_check_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (length(CRAN_check_limit) > 0 && CRAN_check_limit == TRUE) cores <- 2
    ## CRAN_check_limit == TRUE because it is not be a boolean sometimes
    else cores <- detectCores()
  }
  if(verbose)
    message('Parallelizing over ', cores, ' cores.\n')
  
  # We may only create a log-file in debug mode, not for the user
  cl <- if(debug) makeCluster(cores, outfile = 'log.txt') else makeCluster(cores)
  
  ## for showing a progressbar we apparently need to use a SNOW cluster
  hasSNOW <- "doSNOW" %in% installed
  if (verbose) {
    if (hasSNOW) {
      # Remark: SNOW is flagged as superceded by CRAN, but there is no
      # current viable alternative to make a progress bar.
      # The following is thus completely optional in case the user
      # happens to have "doSNOW" installed.
      # We wait for Henrik Bengtsson's work on the "progressr"-API
      do.call("require", list("doSNOW"))
      do.call("registerDoSNOW", list(cl))
      message("Optimization Progress")
      pb <- txtProgressBar(max = length(relevant_families), style = 3)
      progress_fn <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress_fn)
    } else {
      message(length(relevant_families), " families are being searched through.",
              if (length(relevant_families) > cores * 3)
              " This can potentially last several minutes.\n
	      Install the package 'doSNOW' for a progress bar.")
      registerDoParallel(cl)
      opts <- c()
    } 
  } else {
    registerDoParallel(cl)
    opts <- c()
  }
  
  i <- NULL ## BNZ: to prevent an issue, seems to be related to parallel. 
            ##      Do not delete!

  output_liste <- foreach(i = 1:length(relevant_families), .packages = c(), 
                          .errorhandling = 'remove', .verbose = FALSE, 
                          .export = all_funs, .inorder = FALSE,
                          .options.snow = opts) %dopar% {
                            
    fam <- relevant_families[[i]]
    t <- Sys.time()
    
   if(debug) {
     message("Current Family: ", fam$family , " from Package: ",
             fam$package)
   }
    
    result_optim <- eval_with_timeout(
      optimParamsDiscrete(data = data,
                          family = fam[c('package', 'family')],
                          family_info = fam$family_info,
                          method = 'MLE', prior = NULL, 
                          log = fam$family_info$log,
                          optim_method = 'L-BFGS-B', 
                          n_starting_points = 1,
                          debug_error = debug, 
                          show_optim_progress = FALSE,
                          on_error_use_best_result = TRUE, 
                          max_discrete_steps = 100, plot = FALSE,
                          discrete_fast = TRUE, timeout = timeout),
      return_value_on_timeout = "TIMEOUT",
      timeout = 1.5*timeout)
    
    if (length(result_optim) == 1 && result_optim == "TIMEOUT") {
      if (debug) message("Timeout occured for Family ", fam$family)
      result_optim <- NULL
    }
    
    if(!is.null(result_optim) && !is.na(result_optim$value) && 
       !is.infinite(result_optim$value)) {
      output <- new('optimParams', family = fam$family,
                    package = fam$package,
                    estimatedValues = result_optim$par,
                    log_lik = result_optim$value,
                    AIC = result_optim$AIC,
                    BIC = result_optim$BIC,
                    AICc = result_optim$AICc) 
      # aim: check whether solution has good loglik 
      # but does not fit nonetheless
      if(do_sanity) {
        sanity_check <- fitting_sanity_check(output, data, 
                                             continuity = continuity, 
                                             sensitivity = sanity)
        output@sanity <- sanity_check
      }
    } else {
      # a bit redundant, but to keep code modular
      output <- new('optimParams', family = fam$family,
		    package = fam$package,
		    estimatedValues = NA_integer_,
		    log_lik = NA_integer_,
		    AIC = NA_integer_,
		    BIC = NA_integer_,
		    AICc = NA_integer_)
      if(do_sanity) {
        sanity_check <- list(hist_check = NA, int_check = NA, L1_check = NA,
			     good = FALSE)
        output@sanity <- sanity_check
      }
    }
    
   if (debug) {
     cat("Elapsed time for family", fam$family, "from package", fam$package,
	 ":", difftime(Sys.time(), t, units = "secs"), "secs\n")
   }
   
    return(output)
  } # end %dopar%
  stopCluster(cl)

  if(do_sanity) {
    # drop if within-sanity-check yields fail
    families <- sapply(output_liste, function(x) x@family)
    keep_within <- sapply(output_liste, function(x) x@sanity$good)
    if(verbose && sum(keep_within) < length(output_liste)) {
      message("\nSanity Check: Unplausible distribution/fit. Dropping families: ",
      paste(families[!keep_within], collapse = ", "))
    }
    output_liste <- output_liste[keep_within]

    # compare L1-distances with each other, drop if outlier
    families <- sapply(output_liste, function(x) x@family)
    L1 <- sapply(output_liste, function(x) x@sanity$L1_check)
    boxplot <- boxplot(L1, plot = FALSE)
    keep_L1 <- !(L1 %in% boxplot$out & L1 > median(L1))
    # median condition to ensure only exceptionally bad fits are filtered out
    if(verbose && sum(keep_L1) < length(output_liste)) {
      message("Sanity Check: Comparatively bad fit. Dropping families: ",
      paste(families[!keep_L1], collapse = ", "))
    }
    output_liste <- output_liste[keep_L1]
  }
  
  r <- new('globalfit', 
           call = deparse(match.call()),
           data = data, 
           continuity = continuity,
           method = method,
           fits = output_liste)
  r <- sort(r, ic = ic)
  
  r@fits <- r@fits[!is.na(sapply(r@fits, function(x) x %@% ic))]
  
  return(r)
}
