###############################################################################
# Authors:
#   Kim-Anh Le Cao,
#   Benoit Gautier,
#
# created: 2015
# last modified: 12-04-2016
#
# Copyright (C) 2015
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
################################################################################


#' mixOmcis Colors
#' 
#' Default mixOmics color pallete.
#' @param num.vector  A numeric vector of mixOmics colors to return
#'
#' @rdname colors

color.mixo <- function(num.vector){
  
  if (is.factor(num.vector))
    num.vector = as.numeric(num.vector)
  
  if (!is.numeric(num.vector)) {
    stop(paste("num.vector has to be numeric", call. = FALSE))
  }
  
  # these are the colors in the logo (the first 3)
  mixo.gray = gray.colors(1, start = 0.76, gamma = 1)
  
  mixo.col = c('#388ECC', # mixOmics logo blue
               '#F68B33', # mixOmics logo orange
               mixo.gray, # mixOmics logo grey
               '#009E73', # shiny dark green
               '#CC79A7', # shiny purple/pink
               '#F0E442', #shiny yellow
               'black',
               '#D55E00', #shiny dark orange
               '#0072B2', #shiny dark blue
               '#999999'  # shiny grey
  )
  
  #-- checking general input parameters --------------------------------------#
  #---------------------------------------------------------------------------#
  
  n = length(num.vector)
  #-- n: check that there are more colors available than requested
  if (isTRUE(num.vector) > length(mixo.col)) {
    stop(paste("We only have a few mix.colors available, n <= ",
               length(mixo.col)),
         call. = FALSE)
  }
  
  if (isTRUE(!is.finite((num.vector))) ||  (n < 1)) {
    stop("'num.vector' must be an integer vector with positive values.",
         call. = FALSE)
  }
  #-- end checking --#
  #------------------#
  
  return(mixo.col[num.vector])
}


