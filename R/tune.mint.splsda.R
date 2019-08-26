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








#' Estimate the parameters of mint.splsda method
#'
#' Computes Leave-One-Group-Out-Cross-Validation (LOGOCV) scores on a
#' user-input grid to determine optimal values for the sparsity parameters in
#' \code{mint.splsda}.
#'
#' This function performs a Leave-One-Group-Out-Cross-Validation (LOGOCV),
#' where each of \code{study} is left out once. It returns a list of variables
#' of \code{X} that were selected on each of the \code{ncomp} components. Then,
#' a \code{\link{mint.splsda}} can be performed with \code{keepX} set as the
#' output \code{choice.keepX}.
#'
#' All component \eqn{1:\code{ncomp}} are tuned, except the first ones for
#' which a \code{already.tested.X} is provided. See examples below.
#'
#' The function outputs the optimal number of components that achieve the best
#' performance based on the overall error rate or BER. The assessment is
#' data-driven and similar to the process detailed in (Rohart et al., 2016),
#' where one-sided t-tests assess whether there is a gain in performance when
#' adding a component to the model. Our experience has shown that in most case,
#' the optimal number of components is the number of categories in \code{Y} -
#' 1, but it is worth tuning a few extra components to check (see our website
#' and case studies for more details).
#'
#' BER is appropriate in case of an unbalanced number of samples per class as
#' it calculates the average proportion of wrongly classified samples in each
#' class, weighted by the number of samples in each class. BER is less biased
#' towards majority classes during the performance assessment.
#'
#' More details about the prediction distances in \code{?predict} and the
#' supplemental material of the mixOmics article (Rohart et al. 2017).
#'
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y Outcome. Numeric vector or matrix of responses (for multi-response
#' models)
#' @param ncomp Number of components to include in the model (see Details).
#' Default to 1
#' @param study grouping factor indicating which samples are from the same
#' study
#' @param test.keepX numeric vector for the different number of variables to
#' test from the \eqn{X} data set
#' @param already.tested.X if \code{ncomp > 1} Numeric vector indicating the
#' number of variables to select from the \eqn{X} data set on the firsts
#' components
#' @param dist only applies to an object inheriting from \code{"plsda"} or
#' \code{"splsda"} to evaluate the classification performance of the model.
#' Should be a subset of \code{"max.dist"}, \code{"centroids.dist"},
#' \code{"mahalanobis.dist"}. Default is \code{"all"}. See
#' \code{\link{predict}}.
#' @param measure Two misclassification measure are available: overall
#' misclassification error \code{overall} or the Balanced Error Rate \code{BER}
#' @param auc if \code{TRUE} calculate the Area Under the Curve (AUC)
#' performance of the model.
#' @param progressBar by default set to \code{TRUE} to output the progress bar
#' of the computation.
#' @param scale boleean. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default value is FALSE
#' @param light.output if set to FALSE, the prediction/classification of each
#' sample for each of \code{test.keepX} and each comp is returned.
#' @return The returned value is a list with components:
#' \item{error.rate}{returns the prediction error for each \code{test.keepX} on
#' each component, averaged across all repeats and subsampling folds. Standard
#' deviation is also output. All error rates are also available as a list.}
#' \item{choice.keepX}{returns the number of variables selected (optimal keepX)
#' on each component.} \item{choice.ncomp}{returns the optimal number of
#' components for the model fitted with \code{$choice.keepX} }
#' \item{error.rate.class}{returns the error rate for each level of \code{Y}
#' and for each component computed with the optimal keepX}
#'
#' \item{predict}{Prediction values for each sample, each \code{test.keepX} and
#' each comp.} \item{class}{Predicted class for each sample, each
#' \code{test.keepX} and each comp.}
#' @author Florian Rohart
#' @seealso \code{\link{mint.splsda}} and http://www.mixOmics.org for more
#' details.
#' @references Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017).
#' MINT: A multivariate integrative approach to identify a reproducible
#' biomarker signature across multiple experiments and platforms. BMC
#' Bioinformatics 18:128.
#'
#' mixOmics article:
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords multivariate dplot
#' @examples
#'
#' data = stemcells$gene
#' type.id = stemcells$celltype
#' exp = stemcells$study
#'
#' res = mint.splsda(X=data,Y=type.id,ncomp=3,keepX=c(10,5,15),study=exp)
#' out = tune.mint.splsda(X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,
#' study=exp,test.keepX=seq(1,10,1))
#'
#' out$choice.ncomp
#' out$choice.keepX
#'
#' \dontrun{
#'
#' out = tune.mint.splsda(X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,
#' study=exp,test.keepX=seq(1,10,1))
#' out$choice.keepX
#'
#' ## only tune component 2 and keeping 10 genes on comp1
#' out = tune.mint.splsda(X=data,Y=type.id,ncomp=2, study=exp,
#' already.tested.X = c(10),
#' test.keepX=seq(1,10,1))
#' out$choice.keepX
#'
#' }
#'
#' @export tune.mint.splsda
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
    if(is_null(X))
    stop("'X'is missing", call. = FALSE)

    X = as.matrix(X)

    if (length(dim(X)) != 2 || !is.numeric(X))
    stop("'X' must be a numeric matrix.", call. = FALSE)


    # Testing the input Y
    if(is_null(Y))
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
    if (is_null(already.tested.X))
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
    if (is_null(study))
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
    dimnames = list(c(test.keepX), c(paste0('comp', comp.real))))
    mat.mean.error = matrix(nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste0('comp', comp.real))))

    error.per.class.mean = matrix(nrow = nlevels(Y), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(levels(Y)), c(paste0('comp', comp.real))))
    error.per.class.sd = matrix(0,nrow = nlevels(Y), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(levels(Y)), c(paste0('comp', comp.real))))


    error.per.study.keepX.opt = matrix(nrow = nlevels(study), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(levels(study)), c(paste0('comp', comp.real))))

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
    names(error.per.class.keepX.opt) = c(paste0('comp', comp.real))
    names(already.tested.X) = c(paste0('comp', 1:ncomp))

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
        names(auc.mean) = c(paste0('comp', comp.real))
        result$auc = auc.mean
    }

    if(light.output == FALSE)
    {
        names(class.all) = names(prediction.all) = c(paste0('comp', comp.real))
        result$predict = prediction.all
        result$class = class.all
    }
    result$measure = measure
    result$call = match.call()

    class(result) = c("tune.mint.splsda","tune.splsda")

    return(result)
}
