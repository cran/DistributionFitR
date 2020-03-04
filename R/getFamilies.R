## Authors 
## Tim Glockner, tiglockn@mail.uni-mannheim.com
## Adrian Heppeler, aheppele@mail.uni-mannheim.de
## Borui Niklas Zhu, bzhu@mail.uni-mannheim.de
##
## Extract distribution families along with infos regarding 
## their parameters from multiple R packages
##
## Copyright (C) 2019 -- 2020 Tim Glockner, Adrian Heppeler, Borui Niklas Zhu
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



iterate_packages <- function(packages) {
## iterate over packages and extract families and params  
  if (length(packages) == 0) return(list())
  
  res <- list()
  
  ## iterate over packages
  for (i in 1:length(packages)) {
    ## package <- packages[i]
    message("Current Package:", packages[i])
    package_content <- getFamily(packages[i])
    message("\tNumber of families:", length(package_content), "\n")
    
    ## no distribution family contained in package
    if (length(package_content) == 0) next
    
    ## iterate over all families within package
    for (j in 1:length(package_content)) {
      
      stopifnot(package_content[[j]]$package == packages[i])
      
      message("Current Family:", package_content[[j]]$family, "\n")
      ## fetch params for each family and add to 
      params <- tryCatch(getParams(package_content[[j]]),
                         error = function(e) {
                           message("Error occured for family ",
                                   package_content[[j]]$family,
                                   "\n", e)
                           NULL
                         })
      if (length(params)!=0)
        res[[length(res) + 1]] <- c(package_content[[j]],
                                    list(family_info = params))
    }
  }
  
  return(res) 
}


construct_package_list <- function(all.packages) {
  
  if (all.packages) {
    all.packages <- as.vector(installed.packages()[,"Package"])
    
  } else {
    ## only select "base" and "recommended" packages
    ins_pck <- installed.packages()
    
    package_filter <- ((ins_pck[,"Priority"] == "base" |
                        ins_pck[,"Priority"] == "recommended")
      & !is.na(ins_pck[,"Priority"]))
    
    all.packages <- as.vector(ins_pck[package_filter,"Package"])
    
    rm(ins_pck)
  }
  
  return(all.packages)
}


write_file <- function(FamilyList, file = "R/all_families.R") {
  for (i in 1:length(FamilyList)) {
    for (j in 1:length(FamilyList[[i]]$family_info)){
      ij = FamilyList[[i]]$family_info[[j]]
      if (is.double(ij)) {
        FamilyList[[i]]$family_info[[j]] <- signif(ij, 5)
      }
      
    }
  }
  dump(list = "FamilyList", file = file)
}

### Case 1: all.packages given: search for distribution families
### Case 1.1 vector of strings (package names) 
###           -> take families from those packages
### Case 1.2 FALSE -> only take recommended / base packages
### Case 1.3 TRUE -> take all installed packages
### Case 2 all.packages missing: Take families saved in the file

getFamilies <- function(all.packages) {

  ## CASE 2: -> Take families as saved in internal variable FamilyList
  if (missing(all.packages)) {
    return(FamilyList)
  }

  ## CASE 1.2 & 1.3
  if (is.logical(all.packages)) {
    all.packages <- construct_package_list(all.packages = all.packages)
  }
  
  ## CASE 1.1 & 1.2 & 1.3
  FamilyList <- iterate_packages(all.packages)
  ##  write_file(FamilyList=FamilyList, file = file)
  return(FamilyList)
}
