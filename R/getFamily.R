## Authors 
## Manuel Hentschel, mahentsc@mail.uni-mannheim.de
##
## Extract distribution families from a single R package
##
## Copyright (C) 2019 -- 2020 Manuel Hentschel
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


getFamily <- function(pkg){
  
  # load package pkg (can potentially lead to errors 
  # if some requirements are not fulfilled)
  load_successful <- !is(try( suppressPackageStartupMessages(
                              library(pkg, character.only = TRUE)), silent = TRUE),
                              "try-error")
  
  if (!load_successful) {
    warning("Error loading package ", pkg, ". Skipping package.")
    return(list())
  }
  
  start_chars <- c("d", "p", "q", "r")
  # first parameters of the d, p, q, r functions
  first_args <- c("x", "q", "p", "n")   
  
  ## all functions starting with r, d, p or q
  possible_dists <- lsf.str(paste0("package:", pkg), 
                            pattern = paste0("^[", 
                                    paste(start_chars, collapse = ""), "]"))
  
  if (length(possible_dists) == 0) return(list())
  
  l <- vector("list", length(start_chars))
  names(l) <- start_chars
  
  # function for checking whether the first argument of fun in first_arg 
  # (used with first_arg = "x", "n",...)
  check_first_param <- function(fun, first_arg) {
    f_arg <- names(formals(fun))[1]
    !(length(f_arg) == 0) && f_arg == first_arg
  }
  
  for (i in 1:length(start_chars)) {
    char <- start_chars[i]
    # all functions starting with char
    subset <- grep(paste0("^", char), possible_dists, value=TRUE)           

    if (length(subset) != 0) {
      # check if all functions have the correct first arg
      valid_idx <- sapply(subset, check_first_param, first_arg = first_args[i])
      # print(valid_idx)
      l[[char]] <- subset[valid_idx]
    }
  }
  
  get_endings <- function(vec) str_sub(vec, start = 2)
  
  l_endings <- lapply(l, get_endings)    # remove the d, p, q, r suffixes
  
  # we definitely need a function for the density starting with d, 
  # as otherwise we cannot evaluate likelihood function
  # so we only take the endings from p, q and r that also appear in d
  for (char in start_chars[-1]) {
    l_endings[[char]] <- intersect(l_endings[[char]], l_endings$d)
  }
  
  # get a frequency table of the endings
  freq <- table(unlist(l_endings))    
  # only take those distributions that have at least 2 functions implemented
  freq <- freq[freq >= 2]                
  
  families <- names(freq)
  
  # list of lists, where each sublist has the form 
  # list(package=some_pkg, family=some_family)
  families <- lapply(families, function(x) list(package = pkg, family = x))
  
  return(families)
}



