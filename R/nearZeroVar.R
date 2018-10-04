#############################################################################################################
# Authors:
#   This function was borrowed from the package caret nzv.R with some enhancements made by
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2014
# last modified: 12-04-2016
#
# Copyright (C) 2014
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
#############################################################################################################


# --------------------------------------
# nearZeroVar
# --------------------------------------

nearZeroVar = function (x, freqCut = 95/5, uniqueCut = 10)
{
    
    if (is.vector(x))
    x = matrix(x, ncol = 1)
    
    freqRatio = apply(x, 2, function(data)
        {
            data = na.omit(data)
        
            if (length(unique(data)) == length(data))
            { # No duplicate
                return(1)
            } else if (length(unique(data)) == 1) { # Same value
                return(0)
            } else {
                t = table(data)
                return(max(t, na.rm = TRUE)/max(t[-which.max(t)], na.rm = TRUE))
            }
        })
    
    lunique = apply(x, 2, function(data) length(unique(data[!is.na(data)])))
    
    percentUnique = 100 * lunique/nrow(x)
    zeroVar = (lunique == 1) | apply(x, 2, function(data) all(is.na(data)))
    
    out = list()
    out$Position = which((freqRatio > freqCut & percentUnique <= uniqueCut) | zeroVar)
    names(out$Position) = NULL
    out$Metrics = data.frame(freqRatio = freqRatio, percentUnique = percentUnique)
    out$Metrics = out$Metrics[out$Position, ]
    return(out)
}
