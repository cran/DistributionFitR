## Authors 
## Moritz Kern, mkern@mail.uni-mannheim.de
## 
## Define S4 Classes
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

setClass(Class = "globalfit",
         slots = c(call ='character',
                   data = "numeric",
                   continuity = 'logical',
                   method = 'character',
                   fits = 'list'))

setClass(Class = "optimParams",
         slots = c(
           family = "character",
           package = "character",
           estimatedValues = "numeric",
           log_lik = "numeric",
           AIC = "numeric",
           BIC = "numeric",
           AICc = "numeric",
           sanity = "list"
         )
)

setClass(Class = "globalfitSummary",
         slots = c(
           call = "character",
           data = "numeric",
           continuity = "logical",
           method = "character",
           fits = "data.frame",
           ic = 'character'
         )
)
