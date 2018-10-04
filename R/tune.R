#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2016
# last modified: 25-08-2016
#
# Copyright (C) 2016
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


# ========================================================================================================
# tun: chose the optimal number of parameters per component on a "method"
# ========================================================================================================


tune = function (method, # choice of "spls", "splsda", "mint.splsda", "rcc", "pca"
X,
Y,
multilevel = NULL,
ncomp,
study, # mint.splsda
test.keepX = c(5, 10, 15), # all but pca, rcc
test.keepY = NULL, # rcc, multilevel
already.tested.X, # all but pca, rcc
already.tested.Y, #multilevel
mode = "regression", # multilevel
nrepeat = 1, #multilevel, splsda
grid1 = seq(0.001, 1, length = 5), # rcc
grid2 = seq(0.001, 1, length = 5), # rcc
validation = "Mfold", # all but pca
folds = 10, # all but pca
dist = "max.dist", # all but pca, rcc
measure = c("BER"), # all but pca, rcc
auc = FALSE,
progressBar = TRUE, # all but pca, rcc
near.zero.var = FALSE, # all but pca, rcc
logratio = "none", # all but pca, rcc
center = TRUE, # pca
scale = TRUE, # mint, splsda
max.iter = 100, #pca
tol = 1e-09, #pca
light.output = TRUE # mint, splsda
)
{
    choice.method = c("spls", "splsda", "mint.splsda", "rcc", "pca")
    method = match.arg(method, choice.method)
    
    if (method == "mint.splsda") {
        message("Calling 'tune.mint.splsda' with Leave-One-Group-Out Cross Validation (nrepeat = 1)")

        if (missing(ncomp))
        ncomp = 1
        
        result = tune.mint.splsda(X = X, Y = Y,
        ncomp = ncomp,
        study = study,
        test.keepX = test.keepX,
        already.tested.X = already.tested.X,
        dist = dist,
        measure = measure,
        auc = auc,
        progressBar = progressBar,
        scale = scale,
        tol = tol,
        max.iter = max.iter,
        near.zero.var = near.zero.var,
        light.output = light.output)
        
    } else if (method == "rcc") {
        message("Calling 'tune.rcc'")
        
        result = tune.rcc(X = X,
        Y = Y,
        grid1 = grid1,
        grid2 = grid2,
        validation = validation,
        folds = folds,
        plot = plot)
        
    } else if (method == "pca") {
        message("Calling 'tune.pca'")

        if (missing(ncomp))
        ncomp = NULL

        result = tune.pca(X = X,
        ncomp = ncomp,
        center = center,
        scale = scale,
        max.iter = max.iter,
        tol = tol)
        
        
    } else if (method == "splsda") {

        message("Calling 'tune.splsda'")

        if (missing(ncomp))
        ncomp = 1
        
        result = tune.splsda (X = X, Y = Y,
        ncomp = ncomp,
        test.keepX = test.keepX,
        already.tested.X = already.tested.X,
        validation = validation,
        folds = folds,
        dist = dist ,
        measure = measure,
        auc = auc,
        progressBar = progressBar,
        max.iter = max.iter,
        near.zero.var = near.zero.var,
        nrepeat = nrepeat,
        logratio = logratio,
        multilevel = multilevel,
        light.output = light.output)
    } else if (method == "spls") {
        if(missing(multilevel))
        {
            stop("Only a multilevel spls can be tuned")
        } else {
            message("Calling 'tune.splslevel' with method = 'spls'")

            if (missing(ncomp))
            ncomp = 1
            if (missing(already.tested.Y))
            already.tested.Y = NULL
            
            result = tune.splslevel(X = X, Y = Y,
            multilevel = multilevel,
            mode = mode,
            ncomp = ncomp, test.keepX = test.keepX, test.keepY = test.keepY,
            already.tested.X = already.tested.X, already.tested.Y = already.tested.Y)
        }
    }
    
    result$call = match.call()
    return(result)
}
