#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#  Leigh Coonan, Queensland Faculty for Advanced Bioinformatics, Australia
#  Fangzhou Yao, Queensland Faculty for Advanced Bioinformatics, Australia
#  Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Part of this script was borrowed from the prcomp function from the Stats package
#
# created: 2009
# last modified:
#
# Copyright (C) 2009
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


pcasvd = function(X,
ncomp = 3,
retx = TRUE,
center = TRUE,
scale = FALSE)
{
    X = as.matrix(X)
    
    # Borrowed from PRCOMP
    result = svd(get("scale", envir = .GlobalEnv)(X, center = center, scale = scale), nu = 0)
    cen = attr(X, "scaled:center")
    sc = attr(X, "scaled:scale")
    if (any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance.")
    
    # If ncomp is not defined by user at the beginning then this algorithm will calculate an ncomp to be used
    if (is.null(ncomp))
    ncomp = min(nrow(X), ncol(X))
    
    # Borrowed from PRCOMP
    if (ncomp < ncol(X))
    result$v = result$v[, 1:ncomp, drop = FALSE]
    
    result$d = result$d/sqrt(max(1, nrow(X) - 1))
    
    if(retx)
    {
        r = list(sdev = result$d, rotation = result$v,
        X = as.matrix(X) %*% result$v,
        center = if (is.null(cen)) FALSE else cen,
        scale = if (is.null(sc)) FALSE else sc)
    }
    
    if(!retx)
    r = list(sdev = result$d, rotation = result$v)
    
    class(r) = c("pca", "prcomp")
    return(invisible(r))
}
