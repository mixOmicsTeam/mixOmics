#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Liquet, Universite de Bordeaux, France.
#
# created: 2011
# last modified: 24-05-2016
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

# ---------------------------------------------
# withinVariation function 
# ---------------------------------------------
withinVariation = function(X, design){ 
	
    # need a matrix for matrix calculations
    X = as.matrix(X)
    rep.measures = factor(design[, 1])
    factors = design[, -1, drop = FALSE] 
    
    if(any(summary(as.factor(rep.measures)) == 1))
    stop("A multilevel analysis can not be performed when at least one some sample is not repeated.")
        
    # calculate the variation
    # ---------------------------
    # added condition for the spls case where the condition is not needed
    # all we need is the rep.measures
    if ((ncol(factors) == 0) | (ncol(factors) == 1))
    {
      message("Splitting the variation for 1 level factor.")
      
      # save sample names for the output
      indiv.names = rownames(X)
      rownames(X) = as.character(rep.measures)
      
      # compute the mean for each unique individual
      # dealing with specific case with only one subject (leave one out case during prediction)
      X.mean.indiv = matrix(apply(X, 2, tapply, rep.measures, mean, na.rm = TRUE), # to deal with only one subject
                           nrow = length(unique(rep.measures)), ncol = dim(X)[2], dimnames = list(levels(as.factor(rep.measures)), colnames(X)))

      # fill the between matrix with those means 
      Xb = X.mean.indiv[as.character(rep.measures), ]
      
      # compute the within matrix as a difference between the original data 
      # and the between matrix
      Xw = X - Xb
      
      # put dimnames
      dimnames(Xw) = list(indiv.names, colnames(X))    
    } else {  # for 2 levels split
      message("Splitting the variation for 2 level factors.")
      
      ###### off set term
      Xm = colMeans(X)
      
      ###### compute the mean within each subject
      Xs = apply(X, 2, tapply, rep.measures, mean, na.rm = TRUE)
      Xs = Xs[rep.measures, ]
      
      # for the first factor
      xbfact1 = apply(X, 2, tapply, paste0(rep.measures, factors[, 1]), mean, na.rm = TRUE)
      xbfact1 = xbfact1[paste0(rep.measures, factors[, 1]), ]
      
      # for the second factor
      xbfact2 = apply(X, 2, tapply, paste0(rep.measures, factors[, 2]), mean, na.rm = TRUE)
      xbfact2 = xbfact2[paste0(rep.measures, factors[, 2]), ]
      
      #### fixed effect
      # for the first factor
      Xmfact1 = apply(X, 2, tapply, factors[, 1], mean, na.rm = TRUE)
      Xmfact1 = Xmfact1[factors[, 1], ]
      
      # for the second factor
      Xmfact2 = apply(X, 2, tapply, factors[, 2], mean, na.rm = TRUE)
      Xmfact2 = Xmfact2[factors[, 2], ]
      
      # compute the within matrix 
      Xw = X + Xs - xbfact1 - xbfact2 + Xmfact1 + Xmfact2
      Xw = sweep(Xw, 2, 2*Xm)  # see formula in refernece Liquet et al.
      
      # put dimnames
      dimnames(Xw) = dimnames(X)  
    }
    
    return(invisible(Xw))
  }

