#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2015
# last modified: 12-07-2016
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
#############################################################################################################



# logratio.transfo
logratio.transfo = function(X,
logratio = "none", # one of ('none','CLR','ILR')
offset = 0)
{
    if (!(logratio %in% c("none", "CLR", "ILR")))
    stop("Choose one of the three following logratio transformation: 'none', 'CLR' or 'ILR'")
 
    if (logratio == 'ILR')
    {
        if (any(class(X) != 'ilr'))
        {   # data are ilr transformed, then the data lose 1 variable, but we'll use V to reconstruct the matrix
            X = ilr.transfo(X, offset = offset)
        }
    }else if (logratio == 'CLR') {
        X = clr.transfo(X, offset = offset)
    }
    #if logratio = "none", do nothing
    
    return(X)
}


# 1 - ilr transform of the data, isoLMR function from robCompositions package, with changes
# -----------------

# KA changed the function to add a min value when many zeroes in data (prob with log and division by 0 otherwise)
ilr.transfo = function(x, fast = TRUE, offset = 0)
{
    if(any(x==0) & offset ==0)
    stop("make sure you use pseudo counts before normalisation to avoid 0 values with log ratio transformation")
    # ilr transformation
    x.ilr = matrix(NA, nrow = nrow(x), ncol = ncol(x)-1)
    D = ncol(x)
    # KA added: a little something to avoid 0 values
    if (fast)
    {
        for (i in 1 : ncol(x.ilr))
        {
            x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(((apply(as.matrix(x[, (i+1) : D, drop = FALSE]), 1, prod) + offset)^(1 / (D-i))) / (x[,i]+ offset))
            #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop = FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
        }
    } else {
        for (i in 1 : ncol(x.ilr))
        {
            x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(apply(as.matrix(x[, (i+1):D]), 1, function(x){exp(log(x))})/(x[, i]+ offset) + offset)
            #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(apply(as.matrix(x[,(i+1):D]), 1, function(x){exp(log(x))})/(x[,i]))
        }
    }
    class(x.ilr) = 'ilr'
    return(as.matrix(x.ilr))
}




# 2 - back transformation from ilr to clr space
# -------------------
clr.backtransfo = function(x)
{
    # construct orthonormal basis
    V = matrix(0, nrow = ncol(x), ncol = ncol(x)-1)
    for( i in 1:ncol(V) )
    {
        V[1:i, i] = 1/i
        V[i+1, i] = (-1)
        V[, i] = V[, i] * sqrt(i/(i+1))
    }
    rownames(V) = colnames(x)
    return(V)
    
}


# CLR transformation
clr.transfo = function(x, offset = 0)
{
    if(any(x==0) & offset ==0)
    stop("make sure you use pseudo counts before normalisation to avoid 0 values with log ratio transformation")

    # KA added
    #offset = min(x[which(x != 0)])*0.01
    
    
    #if (dim(x)[2] < 2) stop("data must be of dimension greater equal 2")
    if (dim(x)[2] == 1)
    {
        res = list(x.clr = x, gm = rep(1, dim(x)[1]))
    } else{
        geometricmean = function (x) {
            #       if (any(na.omit(x == 0)))
            #         0
            #       else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
            #     }
            # KA changed to
            exp(mean(log(x + offset)))
        }
        gm = apply(x, 1, geometricmean)
        # KA changed
        x.clr = log((x + offset) / (gm))
        res = x.clr #list(x.clr = x.clr, gm = gm)
    }
    class(res) = "clr"
    return(res)
}

