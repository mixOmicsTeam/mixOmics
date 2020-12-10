#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#
# created: 2013
# last modified: 05-10-2017
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


# ========================================================================================================
# tune.spls: chose the optimal number of parameters per component on a spls method, based on MSE
# ========================================================================================================

# I start with no selection on Y. Otherwise we need to be able to tell which submodel of Y is better. I'm afraid a sum(MSE_Yi) is gonna lead to very sparse model (only 1Y), to test.

#' Tuning functions for sPLS method
#' 
#' Computes M-fold or Leave-One-Out Cross-Validation scores on a user-input
#' grid to determine optimal values for the sparsity parameters in \code{spls}.
#' 
#' 
#' This tuning function should be used to tune the parameters in the
#' \code{spls} function (number of components and the number of variables in
#' \code{keepX} to select).
#' 
#' If \code{validation = "loo"}, leave-one-out cross-validation is performed.
#' By default \code{folds} is set to the number of unique individuals. If
#' \code{validation = "Mfold"}, M-fold cross-validation is performed. How many
#' folds to generate is selected by specifying the number of folds in
#' \code{folds}.
#' 
#' Four measures of accuracy are available: Mean Absolute Error (\code{MAE}),
#' Mean Square Error(\code{MSE}), \code{Bias} and \code{R2}. Both MAE and MSE
#' average the model prediction error. MAE measures the average magnitude of
#' the errors without considering their direction. It is the average over the
#' fold test samples of the absolute differences between the Y predictions and
#' the actual Y observations. The MSE also measures the average magnitude of
#' the error. Since the errors are squared before they are averaged, the MSE
#' tends to give a relatively high weight to large errors. The Bias is the
#' average of the differences between the Y predictions and the actual Y
#' observations and the R2 is the correlation between the predictions and the
#' observations. All those measures are averaged across all Y variables in the
#' PLS2 case. We are still improving the function to tune an sPLS2 model,
#' contact us for more details and examples.
#' 
#' The function outputs the optimal number of components that achieve the best
#' performance based on the chosen measure of accuracy. The assessment is
#' data-driven and similar to the process detailed in (Rohart et al., 2016),
#' where one-sided t-tests assess whether there is a gain in performance when
#' adding a component to the model.
#' 
#' See also \code{?perf} for more details.
#' 
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y \code{if(method = 'spls')} numeric vector or matrix of continuous
#' responses (for multi-response models) \code{NA}s are allowed.
#' @param ncomp the number of components to include in the model.
#' @param test.keepX numeric vector for the different number of variables to
#' test from the \eqn{X} data set
#' @param already.tested.X Optional, if \code{ncomp > 1} A numeric vector
#' indicating the number of variables to select from the \eqn{X} data set on
#' the firsts components.
#' @param validation character.  What kind of (internal) validation to use,
#' matching one of \code{"Mfold"} or \code{"loo"} (see below). Default is
#' \code{"Mfold"}.
#' @param folds the folds in the Mfold cross-validation. See Details.
#' @param measure One of \code{MSE} (Mean Squared Error), 
#' \code{MAE} (Mean Absolute Error: MSE without the square), 
#' \code{Bias} (average of the differences), 
#' \code{MAPE} (average of the absolute errors,
#'  as a percentage of the actual values) or \code{R2}. 
#'  Default to \code{MSE}. See details.
#' @param scale Logical. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param progressBar by default set to \code{TRUE} to output the progress bar
#' of the computation.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var Logical, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default value is FALSE
#' @param nrepeat Number of times the Cross-Validation process is repeated.
#' @param multilevel Design matrix for multilevel analysis (for repeated
#' measurements) that indicates the repeated measures on each individual, i.e.
#' the individuals ID. See Details.
#' @param light.output if set to FALSE, the prediction/classification of each
#' sample for each of \code{test.keepX} and each comp is returned.
#' @param cpus Number of cpus to use. If greater than 1, the code is run in
#' parallel.
#' @return A list that contains: \item{error.rate}{returns the prediction error
#' for each \code{test.keepX} on each component, averaged across all repeats
#' and subsampling folds. Standard deviation is also output. All error rates
#' are also available as a list.} \item{choice.keepX}{returns the number of
#' variables selected (optimal keepX) on each component.}
#' \item{choice.ncomp}{returns the optimal number of components for the model
#' fitted with \code{$choice.keepX} and \code{$choice.keepY} }
#' \item{measure}{reminds which criterion was used} \item{predict}{Prediction
#' values for each sample, each \code{test.keepX,test.keepY}, each comp and
#' each repeat. Only if light.output=FALSE}
#' @author Kim-Anh Lê Cao, Benoit Gautier, Francois Bartolo, Florian Rohart,
#' Al J Abadi
#' @seealso \code{\link{splsda}}, \code{\link{predict.splsda}} and
#' http://www.mixOmics.org for more details.
#' @references mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' 
#' PLS and PLS citeria for PLS regression: Tenenhaus, M. (1998). La regression
#' PLS: theorie et pratique. Paris: Editions Technic.
#' 
#' Chavent, Marie and Patouille, Brigitte (2003). Calcul des coefficients de
#' regression et du PRESS en regression PLS1. Modulad n, 30 1-11. (this is the
#' formula we use to calculate the Q2 in perf.pls and perf.spls)
#' 
#' Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction
#' (MSEP) Estimates for Principal Component Regression (PCR) and Partial Least
#' Squares Regression (PLSR). Journal of Chemometrics 18(9), 422-429.
#' 
#' sparse PLS regression mode:
#' 
#' Lê Cao, K. A., Rossouw D., Robert-Granie, C. and Besse, P. (2008). A sparse
#' PLS for variable selection when integrating Omics data. Statistical
#' Applications in Genetics and Molecular Biology 7, article 35.
#' 
#' One-sided t-tests (suppl material):
#' 
#' Rohart F, Mason EA, Matigian N, Mosbergen R, Korn O, Chen T, Butcher S,
#' Patel J, Atkinson K, Khosrotehrani K, Fisk NM, Lê Cao K-A&, Wells CA&
#' (2016). A Molecular Classification of Human Mesenchymal Stromal Cells. PeerJ
#' 4:e1845.
#' @keywords regression multivariate
#' @export
#' @examples
#' 
#' data(liver.toxicity)
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#' 
#' \dontrun{
#' tune = tune.spls(X, Y, ncomp=4, test.keepX = c(5,10,15), measure = "MSE",
#' nrepeat=3, progressBar = TRUE)
#' 
#' tune$choice.ncomp
#' tune$choice.keepX
#' 
#' # plot the results
#' plot(tune)
#' }
tune.spls <- 
    function(X,
             Y,
             test.keepX = NULL,
             test.keepY = NULL,
             ncomp,
             validation = c('Mfold', 'loo'),
             nrepeat = 1,
             folds,
             mode = c('regression', 'canonical', 'classic'),
             measure.tune = c('cor', 'RSS'), ## only if spls model
             BPPARAM = SerialParam(),
             progressBar = FALSE
             ) {
        out = list()
        mode <- match.arg(mode)
        
        X <- .check_numeric_matrix(X, block_name = 'X')
        Y <- .check_numeric_matrix(Y, block_name = 'Y')
        check_cv <- .check_cv_args(validation = validation, 
                                   nrepeat = nrepeat, folds = folds, 
                                   N = nrow(X))
        validation <- check_cv$validation
        nrepeat <- check_cv$nrepeat
        folds <- check_cv$folds
        
        spls.model <- !is.null(test.keepX) | !is.null(test.keepY)
        
        test.keepX <- .change_if_null(arg = test.keepX, default = ncol(X))
        test.keepY <- .change_if_null(arg = test.keepY, default = ncol(Y))
        
        test.keepX <- unique(test.keepX)
        test.keepY <- unique(test.keepY)
        
        if (spls.model) {
            comps <- seq_len(ncomp)
            # TODO check test.keepX and test.keepY
            measure.tune <- match.arg(measure.tune, choices = c('cor', 'RSS'))
            cor.tpred <- 
                cor.upred <- 
                RSS.tpred <- 
                RSS.upred <- 
                array(dim      =    c(length(test.keepX), 
                                      length(test.keepY), 
                                      nrepeat),
                      dimnames = list(paste0('keepX_', test.keepX),
                                      paste0('keepY_', test.keepY),
                                      paste0('repeat_', seq_len(nrepeat))))
        } else {
            if ((test.keepX != ncol(X)) | (test.keepY != ncol(Y)))
                stop("'test.keepX' and 'test.keepY' can only be provided with method = 'spls'", call. = FALSE)
            comps <- ncomp
            test.keepX <- ncol(X)
            test.keepY <- ncol(Y)
            cor.tpred = cor.upred = RSS.tpred = RSS.upred = matrix(nrow = ncomp, ncol = nrepeat,
                                                                   dimnames = list(paste0('comp_', seq_len(ncomp)), paste0('repeat_', seq_len(nrepeat))))
            Q2.tot.ave = matrix(nrow = ncomp, ncol = nrepeat,
                                dimnames = list(paste0('comp_', seq_len(ncomp)), paste0('repeat_', seq_len(nrepeat))))

        }
        choice.keepX = choice.keepY = NULL
        cor.pred = RSS.pred = list()
        cov.pred = list()
        .tune.spls.repeat <- function(test.keepX, test.keepY, X, Y, choice.keepX, choice.keepY, comp, mode, validation, folds)
        {
            out <- list()
            for(keepX in 1:length(test.keepX)){
                for(keepY in 1:length(test.keepY)){
                    # sPLS model, updated with the best keepX
                    pls.res = spls(X = X, Y = Y, 
                                   keepX = c(choice.keepX, test.keepX[keepX]), 
                                   keepY = c(choice.keepY, test.keepY[keepY]), 
                                   ncomp = comp, mode = mode)
                    # fold CV
                    res.perf <- .perf.mixo_pls_cv(pls.res, validation = 'Mfold', folds = folds)
                    out[[paste0('keepX_', keepX)]][[paste0('keepY_', keepY)]] <- list(
                        t.pred.cv = res.perf$t.pred.cv,
                        u.pred.cv = res.perf$u.pred.cv,
                        X.variates =  pls.res$variates$X,
                        Y.variates =  pls.res$variates$Y,
                        Q2.total = res.perf$Q2.total
                    )
                }
            }
            out
        }
        use_progressBar <- progressBar & (is(BPPARAM, 'SerialParam'))
        
            for (comp in comps){
                if (use_progressBar) {
                    msg <- if (spls.model) sprintf("\ntuning component: %s\n", comp) else sprintf("\ntuning pls model ...\n")
                    cat(msg)
                }
                
                cv.repeat.res <- bplapply(seq_len(nrepeat), 
                                        FUN = function(k){ 
                                            if (use_progressBar) {
                                                .progressBar(k/nrepeat)
                                            }
                                            .tune.spls.repeat(test.keepX = test.keepX, test.keepY = test.keepY, X = X, Y = Y, 
                                                              choice.keepX = choice.keepX, choice.keepY = choice.keepY, comp = comp, 
                                                              mode = mode, validation = 'Mfold', folds = folds)}, BPPARAM = BPPARAM)
   
                for(k in seq_len(nrepeat)){
                    for(keepX in 1:length(test.keepX)){
                        for(keepY in 1:length(test.keepY)){
                            t.pred.cv <-  cv.repeat.res[[k]][[paste0('keepX_', keepX)]][[paste0('keepY_', keepY)]]$t.pred.cv
                            u.pred.cv <-  cv.repeat.res[[k]][[paste0('keepX_', keepX)]][[paste0('keepY_', keepY)]]$u.pred.cv
                            X.variates <-  cv.repeat.res[[k]][[paste0('keepX_', keepX)]][[paste0('keepY_', keepY)]]$ X.variates
                            Y.variates <-  cv.repeat.res[[k]][[paste0('keepX_', keepX)]][[paste0('keepY_', keepY)]]$Y.variates
                            Q2.total <- cv.repeat.res[[k]][[paste0('keepX_', keepX)]][[paste0('keepY_', keepY)]]$Q2.total
                            
                            if (spls.model)
                            {
                                # extract the predicted components: 
                                # if(measure.tune == 'cor' ){
                                    cor.tpred[keepX, keepY, k] = abs(cor(t.pred.cv[, comp], X.variates[, comp]))
                                    cor.upred[keepX, keepY,k] = abs(cor(u.pred.cv[, comp], Y.variates[, comp]))
                                # }
                                # if(measure.tune == 'RSS'){
                                    # RSS: no abs values here
                                    RSS.tpred[keepX, keepY, k] = sum((t.pred.cv[, comp] - X.variates[, comp])^2)/(nrow(X) -1) 
                                    RSS.upred[keepX, keepY, k] = sum((u.pred.cv[, comp] - Y.variates[, comp])^2)/(nrow(X) -1)
                                # }
                                # covariance between predicted variates
                                ##cov.variate.pred[keepX, keepY, k] = cov(t.pred.cv[, comp], u.pred.cv[, comp])
                                
                            } else {
                                
                                # extract Q2.total for a PLS, we could extract other outputs such as R2, MSEP etc (only valid for regression)
                                Q2.tot.ave[, k] = Q2.total 
                                
                                # extract the predicted components per dimension, take abs value
                                cor.tpred[, k] = diag(abs(cor(t.pred.cv, X.variates)))
                                cor.upred[, k] = diag(abs(cor(u.pred.cv, Y.variates)))
                                
                                # RSS: no abs values here
                                RSS.tpred[, k] = apply((t.pred.cv - X.variates)^2, 2, sum)/(nrow(X) -1)
                                RSS.upred[, k] = apply((u.pred.cv - Y.variates)^2, 2, sum)/(nrow(X) -1)
                            }
    
                        } # end keepY
                    } #end keepX
                } # end repeat
                cat('\t')
                
                ## add mean and sd across repeats for output
                .get_mean_and_sd <- function(arr) {
                    list(values = arr, 
                         mean   = apply(arr, c(1,2), mean), 
                         sd     = apply(arr, c(1,2), sd))
                } 
                cor.pred$u[[paste0('comp_', comp)]] = .get_mean_and_sd(cor.upred)
                cor.pred$t[[paste0('comp_', comp)]] = .get_mean_and_sd(cor.tpred)
                RSS.pred$u[[paste0('comp_', comp)]] = .get_mean_and_sd(RSS.upred)
                RSS.pred$t[[paste0('comp_', comp)]] = .get_mean_and_sd(RSS.tpred)
                
                t.test.arr <- function(arr, is_cor) {
                    if (dim(arr)[3] < 3) ## low nrepeat, no t.test
                    {
                        extremum <- ifelse(is_cor, max, min)
                        ind.opt <- which(arr == extremum(arr), arr.ind = TRUE)
                        
                        return(c(ind.choice.keepX = ind.opt[1], 
                                 ind.choice.keepY = ind.opt[2]))
                        
                    }
                    choice.keepX_i <- 1
                    choice.keepY_j <- 1
                    for (keepX_i in seq_len(dim(arr)[1])[-1]) {
                        for (keepY_j in seq_len(dim(arr)[2])[-1])
                        {
                            x <- arr[choice.keepX_i, choice.keepY_j, ]
                            y <- arr[keepX_i, keepY_j, ]
                            t.test.res <- t.test(x,
                                                 y,
                                                 alternative = ifelse(is_cor, 'less', 'greater'),
                                                 paired = FALSE)
                            if (t.test.res$p.value < 0.05)
                            {
                                choice.keepX_i <- keepX_i
                                choice.keepY_j <- keepY_j
                            }
                            
                        }
                    }
                    return(c(ind.choice.keepX = choice.keepX_i, 
                             ind.choice.keepY = choice.keepY_j))
                }
                
                
                if (spls.model)
                {
                # choose the best keepX and keepY based on type.tune
                if(mode != 'canonical'){  #regression, invariant, classic
                    # define best keepX and keepY based on u
                    if(measure.tune == 'cor'){
                        cor.component = cor.pred$u[[paste0('comp_', comp)]]
                        # index = which(cor.component == max(cor.component), arr.ind = TRUE)
                        index <- t.test.arr(arr = cor.component$values, is_cor = TRUE)
                    }else{ # if type.tune = 'RSS'
                        RSS.component = RSS.pred$u[[paste0('comp_', comp)]]
                        # index = which(RSS.component == min(RSS.component), arr.ind = TRUE)
                        index <- t.test.arr(arr = RSS.component$values, is_cor = FALSE)
                    }
                    choice.keepX = c(choice.keepX, test.keepX[index['ind.choice.keepX']])
                    choice.keepY = c(choice.keepY, test.keepY[index['ind.choice.keepY']])
                    
                }else{  # mode = 'canonical'
                    if(measure.tune == 'cor'){
                        cor.component.t = cor.pred$t[[paste0('comp_', comp)]]
                        cor.component.u = cor.pred$u[[paste0('comp_', comp)]]
                        index.t = t.test.arr(arr = cor.component.t$values, is_cor = TRUE)
                        index.u = t.test.arr(arr = cor.component.u$values, is_cor = TRUE)
                    }else{ # if type.tune = 'RSS'
                        RSS.component.t = RSS.pred$t[[paste0('comp_', comp)]]
                        RSS.component.u = RSS.pred$u[[paste0('comp_', comp)]]
                        index.t = t.test.arr(arr = RSS.component.t$values, is_cor = FALSE)
                        index.u = t.test.arr(arr = RSS.component.u$values, is_cor = FALSE)
                    }
                    choice.keepX = c(choice.keepX, test.keepX[index.t['ind.choice.keepX']])
                    choice.keepY = c(choice.keepY, test.keepY[index.u['ind.choice.keepY']])
                } # canonical
                }
                    } # end comp
            if (spls.model)
            {
               out$choice.keepX = choice.keepX
               out$choice.keepY = choice.keepY
            }

        # } # end sPLS
        
        if(spls.model){
            out$cor.pred = cor.pred
            out$RSS.pred = RSS.pred
        }else{
            out$cor.pred = cor.pred
            out$RSS.pred = RSS.pred
            out$Q2.tot.ave = apply(Q2.tot.ave, 1, mean)
        }
        ### evaluate all for output except X and Y to save memory
        ## eval all but X and Y
        mc <- mget(names(formals())[-1:-2], sys.frame(sys.nframe()))
        ## replace function, X and Y with unevaluated call
        mc <- as.call(c(as.list(match.call())[1:3], mc))
        out <- c(list(call = mc), out)
        class(out) <- if (spls.model) c('tune.pls', 'tune.spls') else c('tune.pls')
        return(out)
    } 
