#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2016
# last modified: 05-10-2017
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
# tune.mint.splsda: chose the optimal number of parameters per component on a mint.splsda method
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 1.
# study: grouping factor indicating which samples are from the same study
# test.keepX: grid of keepX among which to chose the optimal one
# already.tested.X: a vector giving keepX on the components that were already tuned
# dist: distance to classify samples. see predict
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# progressBar: show progress,
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# nrepeat: number of replication of the Mfold process
# light.output: if TRUE, only the most important outputs are given (and calculated)


tune.mint.splsda = function (X, Y,
ncomp = 1,
study,
test.keepX = c(5, 10, 15),
already.tested.X,
dist = "max.dist",
measure = "BER", # one of c("overall","BER")
auc = FALSE,
progressBar = TRUE,
scale = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
light.output = TRUE # if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
)
{    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #------------------#
    #-- check entries --#
    if(missing(X))
    stop("'X'is missing", call. = FALSE)

    X = as.matrix(X)

    if (length(dim(X)) != 2 || !is.numeric(X))
    stop("'X' must be a numeric matrix.", call. = FALSE)
    
    
    # Testing the input Y
    if(missing(Y))
    stop("'Y'is missing", call. = FALSE)
    if (is.null(Y))
    stop("'Y' has to be something else than NULL.", call. = FALSE)
    
    if (is.null(dim(Y)))
    {
        Y = factor(Y)
    }  else {
        stop("'Y' should be a factor or a class vector.", call. = FALSE)
    }
    
    if (nlevels(Y) == 1)
    stop("'Y' should be a factor with more than one level", call. = FALSE)
    
    
    
    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    
    #-- measure
    if (length(measure) > 1)
    measure = measure[1]
    
    if (!measure %in% c("overall", "BER"))
    stop("'measure must be one of 'overall' or 'BER'", call. = FALSE)
    
    #if ((!is.null(already.tested.X)) && (length(already.tested.X) != (ncomp - 1)) )
    #stop("The number of already tested parameters should be NULL or ", ncomp - 1, " since you set ncomp = ", ncomp)
    if (missing(already.tested.X))
    {
        already.tested.X = NULL
    } else {
            if(is.list(already.tested.X))
            stop("''already.tested.X' must be a vector of keepX values")
            
            message(paste("Number of variables selected on the first", length(already.tested.X), "component(s):", paste(already.tested.X,collapse = " ")))
    }
    if(length(already.tested.X) >= ncomp)
    stop("'ncomp' needs to be higher than the number of components already tuned, which is length(already.tested.X)=",length(already.tested.X) , call. = FALSE)
    

    # -- check using the check of mint.splsda
    Y.mat = unmap(Y)
    colnames(Y.mat) = levels(Y)
    
    check = Check.entry.pls(X, Y = Y.mat, ncomp = ncomp, mode="regression", scale=scale,
    near.zero.var=near.zero.var, max.iter=max.iter ,tol=tol ,logratio="none" ,DA=TRUE, multilevel=NULL)
    X = check$X
    ncomp = check$ncomp
    
    
    # -- study
    #set the default study factor
    if (missing(study))
    stop("'study' is missing", call. = FALSE)

    if (length(study) != nrow(X))
    stop(paste0("'study' must be a factor of length ",nrow(X),"."))
    
    if (any(table(study) <= 1))
    stop("At least one study has only one sample, please consider removing before calling the function again", call. = FALSE)
    if (any(table(study) < 5))
    warning("At least one study has less than 5 samples, mean centering might not do as expected")
    
    if(sum(apply(table(Y,study)!=0,2,sum)==1) >0)
    stop("At least one study only contains a single level of the multi-levels outcome Y. The MINT algorithm cannot be computed.")
    
    if(sum(apply(table(Y,study)==0,2,sum)>0) >0)
    warning("At least one study does not contain all the levels of the outcome Y. The MINT algorithm might not perform as expected.")


    #-- dist
    
    choices = c("max.dist", "centroids.dist", "mahalanobis.dist")
    dist = match.arg(dist, choices, several.ok = FALSE)
    
    #-- light.output
    if (!is.logical(light.output))
    stop("'light.output' must be either TRUE or FALSE", call. = FALSE)
    
    #-- test.keepX
    if (is.null(test.keepX) | length(test.keepX) == 1 | !is.numeric(test.keepX))
    stop("'test.keepX' must be a numeric vector with more than two entries", call. = FALSE)

    #-- end checking --#
    #------------------#
    
    
    
    #-- cross-validation approach  ---------------------------------------------#
    #---------------------------------------------------------------------------#
    test.keepX = sort(test.keepX) #sort test.keepX so as to be sure to chose the smallest in case of several minimum
    
    # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
    if ((!is.null(already.tested.X)))
    {
        comp.real = (length(already.tested.X) + 1):ncomp
    } else {
        comp.real = 1:ncomp
    }

    mat.error = matrix(nrow = length(test.keepX), ncol = 1,
    dimnames = list(test.keepX,1))
    rownames(mat.error) = test.keepX
    
    error.per.class = list()
    
    mat.sd.error = matrix(0,nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real))))
    mat.mean.error = matrix(nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real))))
    
    error.per.class.mean = matrix(nrow = nlevels(Y), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(levels(Y)), c(paste('comp', comp.real))))
    error.per.class.sd = matrix(0,nrow = nlevels(Y), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(levels(Y)), c(paste('comp', comp.real))))
    
    
    error.per.study.keepX.opt = matrix(nrow = nlevels(study), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(levels(study)), c(paste('comp', comp.real))))
    
    if(light.output == FALSE)
    prediction.all = class.all = list()
    if(auc)
    auc.mean=list()
    
    error.per.class.keepX.opt=list()
    
    # successively tune the components until ncomp: comp1, then comp2, ...
    for(comp in 1:length(comp.real))
    {
        
        if (progressBar == TRUE)
        cat("\ncomp",comp.real[comp], "\n")
        
        result = LOGOCV (X, Y, ncomp = 1 + length(already.tested.X), study = study,
        choice.keepX = already.tested.X,
        test.keepX = test.keepX, measure = measure,
        dist = dist, near.zero.var = near.zero.var, progressBar = progressBar, scale = scale, max.iter = max.iter, auc = auc)
        

        # in the following, there is [[1]] because 'tune' is working with only 1 distance and 'MCVfold.splsda' can work with multiple distances
        mat.mean.error[, comp]=result[[measure]]$error.rate.mean[[1]]
        if (!is.null(result[[measure]]$error.rate.sd[[1]]))
        mat.sd.error[, comp]=result[[measure]]$error.rate.sd[[1]]
        
        # confusion matrix for keepX.opt
        error.per.class.keepX.opt[[comp]]=result[[measure]]$confusion[[1]]
        
        # best keepX
        already.tested.X = c(already.tested.X, result[[measure]]$keepX.opt[[1]])

        # error per study for keepX.opt
        error.per.study.keepX.opt[,comp] = result[[measure]]$error.per.study.keepX.opt[[1]]
        
        if(light.output == FALSE)
        {
            #prediction of each samples for each fold and each repeat, on each comp
            class.all[[comp]] = result$class.comp[[1]]
            prediction.all[[comp]] = result$prediction.comp
        }
        if(auc)
        auc.mean[[comp]]=result$auc
        
    } # end comp
    names(error.per.class.keepX.opt) = c(paste('comp', comp.real))
    names(already.tested.X) = c(paste('comp', 1:ncomp))
    
    if (progressBar == TRUE)
    cat('\n')
    
    # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
    if(nlevels(study) > 2 & length(comp.real) >1)
    {
        opt = t.test.process(error.per.study.keepX.opt)
        ncomp_opt = comp.real[opt]
    } else {
        ncomp_opt = NULL
    }
    
    result = list(
    error.rate = mat.mean.error,
    choice.keepX = already.tested.X,
    choice.ncomp = list(ncomp = ncomp_opt, values = error.per.study.keepX.opt),
    error.rate.class = error.per.class.keepX.opt)
    
    if(auc)
    {
        names(auc.mean) = c(paste('comp', comp.real))
        result$auc = auc.mean
    }
    
    if(light.output == FALSE)
    {
        names(class.all) = names(prediction.all) = c(paste('comp', comp.real))
        result$predict = prediction.all
        result$class = class.all
    }
    result$measure = measure
    result$call = match.call()

    class(result) = c("tune.mint.splsda","tune.splsda")

    return(result)
}
