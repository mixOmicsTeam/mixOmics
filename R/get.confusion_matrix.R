################################################################################
# Authors:
#   Florian Rohart,
#
# created: 2017
# last modified: 31-03-2017
#
# Copyright (C) 2017
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


# =============================================================================
# get.confusion_matrix: create confusion table between a vector of true classes
#   and a vector of predicted classes
# =============================================================================

# truth: factor of true classes
# levels: levels of the 'truth' factor. Optional parameters that can be used
#   when there are some missing levels in `truth' compared to the fitted model
# predicted: vector of predicted classes. Can contain NA.

get.confusion_matrix = function(truth, all.levels, predicted)
{
    
    if(length(truth) != length(predicted))
    stop("'truth' and 'predicted' must be of same length")
    
    if(!is.factor(truth))
    truth = factor(truth)
    
    if(missing(all.levels))
    all.levels = levels(truth)
    
    #print(all.levels)
    
    nlevels.truth = length(all.levels)
    
    ClassifResult = array(0,c(nlevels.truth, nlevels.truth + 1))
    #+1 for NA prediction
    rownames(ClassifResult) = all.levels
    colnames(ClassifResult) = paste("predicted.as.",c(all.levels, "NA"),
    sep = "")
    
    #--------record of the classification accuracy for each level of Y
    for(i in 1:nlevels.truth)
    {
        ind.i = which(truth == all.levels[i])
        for(ij in 1:nlevels.truth)
        ClassifResult[i,ij] = sum(predicted[ind.i] == all.levels[ij],
        na.rm = TRUE)
        
        # if some NA, we add them in the last column (ij+1 = nlevels.truth + 1)
        if(sum(is.na(predicted[ind.i]))>0)
        ClassifResult[i,ij+1] =  sum(is.na(predicted[ind.i]))
    }
    
    # if no NA in the prediction, we remove the last column
    if(sum(is.na(predicted))==0)
    ClassifResult = ClassifResult [, -(nlevels.truth+1)]

    ClassifResult
}

# calculate BER from a confusion matrix
get.BER = function(confusion)
{
    #if(!is.numeric(X)| !is.matrix(X) | length(dim(X)) != 2 | nrow(X)!=ncol(X))
    #stop("'X' must be a square numeric matrix")
    
    nlev = nrow(confusion)
    #calculation of the BER
    ClassifResult.temp = confusion
    diag(ClassifResult.temp) = 0
    BER = sum(apply(ClassifResult.temp,1,sum,na.rm = TRUE)/apply(confusion,1,
        sum,na.rm = TRUE),na.rm = TRUE)/nlev
    return(BER)
}

