#############################################################################################################
# Author :
#   Kim-Anh Le Cao, ARC Centre of Excellence in Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Leigh Coonan, Student, University of Quuensland, Australia
#   Fangzhou Yao, Student, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2011
# last modified: 21-04-2016
#
# Copyright (C) 2011
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

tune.pca =
function(X,
ncomp = NULL,
center = TRUE, 	# sets the mean of the data to zero, ensures that the first PC describes the direction of the maximum variance
scale = FALSE, 	# variance is unit accross different units
max.iter = 500,
tol = 1e-09,
logratio = 'none',# one of ('none','CLR','ILR')
V = NULL,
multilevel = NULL)
{
    
    
    result = pca(X = X, ncomp = ncomp,
    center = center, scale = scale,
    max.iter = max.iter, tol = tol,
    logratio = logratio, V = V,
    multilevel = multilevel)
    
    is.na.X = is.na(X)
    na.X = FALSE
    if (any(is.na.X)) na.X = TRUE
    
    #  list eigenvalues, prop. of explained varience and cumulative proportion of explained variance
    prop.var = result$explained_variance
    cum.var = result$cum.var
    
    ind.show = min(10, ncomp)
    
    print(result)

    # Plot the principal components and explained variance
    # note: if NA values, we have an estimation of the variance using NIPALS
    if(!na.X)
    {
        ylab = "Proportion of Explained Variance"
    } else{
        ylab = "Estimated Proportion of Explained Variance"
    }
    barplot(prop.var[1:result$ncomp], names.arg = 1:result$ncomp, xlab = "Principal Components",
    ylab = ylab)
    
    result$call = match.call()
    
    class(result) = "tune.pca"
    return(invisible(result))
}
