################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#
# created: 2015
# last modified:
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


#-- spectral colors --#
#---------------------#
#' @importFrom grDevices colorRamp
color.spectral =
  function (n, alpha = 1)
  {
    #-- checking general input parameters -------------------------------------#
    #--------------------------------------------------------------------------#

    #-- n
    if (length(n) > 1 || !is.finite(n))
    stop("'n' must be an integer positive value.", call. = FALSE)

    if (n < 1)
    stop("'n' must be an integer positive value.", call. = FALSE)

    #-- alpha
    if (length(alpha) > 1 || !is.finite(alpha))
    stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)

    if (alpha < 0 || alpha > 1)
    stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)

    alpha = round(255 * alpha)

    #-- end checking --#
    #------------------#

    ramp = colorRamp(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
                       "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F",
                       "#9E0142"), space = "Lab")

    rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
  }

