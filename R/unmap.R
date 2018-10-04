#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
#
# This function was borrowed from the mclust package and modified for mixOmics
#
# created: 2013
# last modified: 12-04-2016
#
# Copyright (C) 2013
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

# ---------------------------------------------------
# unmap variates.A variable for (s)plsda
# ---------------------------------------------------

unmap = function (classification, groups = NULL, noise = NULL)
{
    n = length(classification)
    u = sort(unique(classification))
    levels =  levels(classification)### Add levels
    
    if (is.null(groups))
    {
        groups = u
    } else {
        if (any(match(u, groups, nomatch = 0) == 0))
        stop("groups incompatible with classification")
        miss = match(groups, u, nomatch = 0) == 0
    }
    
    cgroups = as.character(groups)
    if (!is.null(noise))
    {
        noiz = match(noise, groups, nomatch = 0)
        if (any(noiz == 0))
        stop("noise incompatible with classification")
        
        groups = c(groups[groups != noise], groups[groups == noise])
        noise = as.numeric(factor(as.character(noise), levels = unique(groups)))
    }
    
    groups = as.numeric(factor(cgroups, levels = unique(cgroups)))
    classification = as.numeric(factor(as.character(classification), levels = unique(cgroups)))
    k = length(groups) - length(noise)
    nam = levels(groups)
    
    if (!is.null(noise))
    {
        k = k + 1
        nam = nam[1:k]
        nam[k] = "noise"
    }
    
    z = matrix(0, n, k, dimnames = c(names(classification), nam))
    for (j in 1:k) z[classification == groups[j], j] = 1
    attr(z, "levels") = levels
    z
}

# ---------------------------------------------------
# map variable for (s)plsda
# ---------------------------------------------------

map = function (Y)
{
    nrowY = nrow(Y)
    cl = numeric(nrowY)
    I = 1:nrowY
    J = 1:ncol(Y)
    for (i in I)
    {
        cl[i] = (J[Y[i, ] == max(Y[i, ])])[1]
    }
    return(cl)
}
