#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Leigh Coonan, Queensland Faculty for Advanced Bioinformatics, Australia
#
# created: 2010
# last modified: 19-04-2016
#
# Copyright (C) 2010
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


plot.pca  = #plot.spca <- plot.ipca <- plot.sipca <-
function(   x,
            ncomp = min(10, length(x$sdev)),
            type = "barplot", # either barplot or any other type available in plot, as "l","b","p",..
            explained.var=TRUE,
            ...)
{
    
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- ncomp
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
    stop("invalid value for 'ncomp'.", call. = FALSE)
    
    ncomp = round(ncomp)
    
    if (ncomp > length(x$sdev))
    stop("'ncomp' must be lower or equal than ", length(x$sdev), ".",
    call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    #-- scree plot -------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    variances = (x$sdev^2)[1:ncomp] # relative variance
    ylab = "Variance"
    if(explained.var==TRUE)
    {
        variances=variances/x$var.tot #explained variances
        ylab = "Explained Variance"
    }
    if (type == "barplot")
    {
        barplot(variances, names.arg = seq(1, ncomp), xlab = "Principal Components", ylab = ylab,...)
    } else {
        plot(variances, type = type, axes = FALSE,
        xlab = "Principal Components",
        ylab = ylab,... )
        axis(1, at = 1:ncomp)
        axis(2)
    }
    
}


