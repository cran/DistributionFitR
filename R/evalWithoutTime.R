## Authors 
## Benedikt Geier, bgeier@mail.uni-mannheim.de
##
## Various help functions
##
## Copyright (C) 2019 -- 2020 Benedikt Geier, Henrik Bengtsson
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



# Function for stopping an expression if it takes too long
# Simplified version of https://github.com/HenrikBengtsson/R.utils/issues/74
eval_with_timeout <- function(expr, envir = parent.frame(), timeout,
                              return_value_on_timeout = NULL) {
  # substitute expression so it is not executed as soon it is used
  expr <- substitute(expr)
  
  # for Borui
  if (.Platform$OS.type != "unix") {
    return(eval(expr, envir = envir))
  }
  
  # execute expr in separate fork
  myfork <- parallel::mcparallel({
    eval(expr, envir = envir)
  }, silent = FALSE)
  
  # wait max n seconds for a result.
  myresult <- parallel::mccollect(myfork, wait = FALSE, timeout = timeout)
  # kill fork after collect has returned
  tools::pskill(myfork$pid, tools::SIGKILL)
  tools::pskill(-1 * myfork$pid, tools::SIGKILL)
  
  # clean up:
  parallel::mccollect(myfork, wait = FALSE)
  
  # timeout?
  if (is.null(myresult)) {
    myresult <- return_value_on_timeout
  }
  
  # move this to distinguish between timeout and NULL returns
  myresult <- myresult[[1]]
  
  if ("try-error" %in% class(myresult)) {
    stop(attr(myresult, "condition"))
  }
  
  # send the buffered response
  return(myresult)
}
