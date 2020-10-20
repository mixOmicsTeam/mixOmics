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
#' @param scale Boolean. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param progressBar by default set to \code{TRUE} to output the progress bar
#' of the computation.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
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
             test.keepX,
             test.keepY,
             ncomp,
             nrepeat,
             mode,
             type.tune, # ! shoud be null for a PLS model
             pls.model) {
        
        
        out = list()
        
        if(isFALSE(pls.model)){
            cor.tpred = cor.upred = RSS.tpred = RSS.upred = array(dim = c(length(test.keepX), length(test.keepY), nrepeat),
                                                                  dimnames = list(paste0('keepX', test.keepX), 
                                                                                  paste0('keepY', test.keepY),
                                                                                  paste0('repeat', 1:nrepeat)))
            best.keepX = best.keepY = NULL
        }else{
            type.tune = NULL
            cor.tpred = cor.upred = RSS.tpred = RSS.upred = matrix(nrow = ncomp, ncol = nrepeat, 
                                                                   dimnames = list(paste0('comp', 1:ncomp), paste0('repeat', 1:nrepeat)))
            Q2.tot.ave = matrix(nrow = ncomp, ncol = nrepeat,
                                dimnames = list(paste0('comp', 1:ncomp), paste0('repeat', 1:nrepeat)))
            
        }
        cor.pred = RSS.pred = list()
        cov.pred = list()
        
        
        
        ## for a PLS only to extract Q2.total (or anything else)
        if(isTRUE(pls.model)){
            for(k in 1:nrepeat){
                cat('repeat', k, '\n')
                
                # run a full PLS model up to the last comp
                pls.res = spls(X = X, Y = Y, 
                               keepX = c(rep(ncol(X), ncomp)), 
                               keepY = c(rep(ncol(Y), ncomp)), 
                               ncomp = ncomp, mode = mode)
                # fold CV
                res.perf = .perf.mixo_pls_folds(pls.res, validation = 'Mfold', folds = 10)
                
                # extract Q2.total for a PLS, we could extract other outputs such as R2, MSEP etc (only valid for regression)
                Q2.tot.ave[, k] = res.perf$Q2.total 
                
                # extract the predicted components per dimension, take abs value
                cor.tpred[, k] = diag(abs(cor(res.perf$t.pred.cv, pls.res$variates$X)))
                cor.upred[, k] = diag(abs(cor(res.perf$u.pred.cv, pls.res$variates$Y)))
                
                # RSS: no abs values here
                RSS.tpred[, k] = apply((res.perf$t.pred.cv - pls.res$variates$X)^2, 2, sum)/(nrow(X) -1)
                RSS.upred[, k] = apply((res.perf$u.pred.cv - pls.res$variates$Y)^2, 2, sum)/(nrow(X) -1)
            } #end repeat       
            
            # # calculate mean across repeats
            cor.pred$u = apply(cor.upred, 1, mean) 
            cor.pred$t = apply(cor.tpred, 1, mean)
            RSS.pred$u = apply(RSS.upred, 1, mean)
            RSS.pred$t = apply(RSS.tpred, 1, mean)
            
        }else{ # if sPLS model 
            for (comp in 1:ncomp){
                cat('Comp', comp, '\n')
                for(k in 1:nrepeat){
                    cat('repeat', k, '\n')
                    for(keepX in 1:length(test.keepX)){
                        #cat('KeepX', list.keepX[keepX], '\n')
                        for(keepY in 1:length(test.keepY)){
                            #cat('KeepY', list.keepY[keepY], '\n')
                            
                            # sPLS model, updated with the best keepX
                            spls.res = spls(X = X, Y = Y, 
                                            keepX = c(best.keepX, test.keepX[keepX]), 
                                            keepY = c(best.keepY, test.keepY[keepY]), 
                                            ncomp = comp, mode = mode)
                            # fold CV
                            res.perf = .perf.mixo_pls_folds(spls.res, validation = 'Mfold', folds = 10)
                            
                            # extract the predicted components: 
                            if(type.tune == 'cor' ){
                                cor.tpred[keepX, keepY, k] = cor(res.perf$t.pred.cv[,comp], spls.res$variates$X[, comp])
                                cor.upred[keepX, keepY,k] = cor(res.perf$u.pred.cv[,comp], spls.res$variates$Y[, comp])
                            }
                            if(type.tune == 'RSS'){
                                # RSS: no abs values here
                                RSS.tpred[keepX, keepY, k] = sum((res.perf$t.pred.cv[,comp] - spls.res$variates$X[, comp])^2)/(nrow(X) -1) 
                                RSS.upred[keepX, keepY, k] = sum((res.perf$u.pred.cv[,comp] - spls.res$variates$Y[, comp])^2)/(nrow(X) -1)
                            }
                            # covariance between predicted variates
                            ##cov.variate.pred[keepX, keepY, k] = cov(res.perf$t.pred.cv[,comp], res.perf$u.pred.cv[,comp])
                        } # end keepY
                    } #end keepX
                } # end repeat
                cat('\t')
                
                # # calculate mean across repeats
                cor.pred$u[[comp]] = apply(cor.upred, c(1,2), mean)  #mean(cor.upred[keepX, keepY,])
                cor.pred$t[[comp]] = apply(cor.tpred, c(1,2), mean)  #mean(cor.tpred[keepX, keepY,])
                RSS.pred$u[[comp]] = apply(RSS.upred, c(1,2), mean)  #mean(RSS.upred[keepX, keepY,])
                RSS.pred$t[[comp]] = apply(RSS.tpred, c(1,2), mean)  #mean(RSS.tpred[keepX, keepY,])
                
                # choose the best keepX and keepY based on type.tune
                if(mode != 'canonical'){  #regression, invariant, classic
                    # define best keepX and keepY based on u
                    if(type.tune == 'cor'){
                        cor.component = cor.pred$u[[comp]]
                        index = which(cor.component == max(cor.component), arr.ind = TRUE)
                    }else{ # if type.tune = 'RSS'
                        RSS.component = RSS.pred$u[[comp]]
                        index = which(RSS.component == min(RSS.component), arr.ind = TRUE)
                    }
                    best.keepX = c(best.keepX, test.keepX[index[1,1]])
                    best.keepY = c(best.keepY, test.keepY[index[1,2]])
                    
                }else{  # mode = 'canonical'
                    if(type.tune == 'cor'){
                        cor.component.t = cor.pred$t[[comp]]
                        cor.component.u = cor.pred$u[[comp]]
                        index.t = which(cor.component.t == max(cor.component.t), arr.ind = TRUE)
                        index.u = which(cor.component.u == max(cor.component.u), arr.ind = TRUE)
                    }else{ # if type.tune = 'RSS'
                        RSS.component.t = RSS.pred$t[[comp]]
                        RSS.component.u = RSS.pred$u[[comp]]
                        index.t = which(RSS.component.t == min(RSS.component.t), arr.ind = TRUE)
                        index.u = which(RSS.component.u == min(RSS.component.u), arr.ind = TRUE)
                    }
                    best.keepX = c(best.keepX, test.keepX[index.t[1,1]])
                    best.keepY = c(best.keepY, test.keepY[index.u[1,2]])
                } # canonical
            } # end comp
            
            out$best.keepX = best.keepX
            out$best.keepY = best.keepY
        } # end sPLS
        
        
        
        
        # 
        # # # summary of results
        # summary = list()
        # for(k in 1:ncomp){
        #   RSS.upred = RSS.pred$u[[k]]
        #   RSS.tpred = RSS.pred$t[[k]]
        # 
        #   cor.upred = cor.pred$u[[k]]
        #   cor.tpred = cor.pred$t[[k]]
        # 
        #   # u.pred
        #   index = which(cor.upred == max(cor.upred), arr.ind = TRUE)
        #   #list.keepX[index[1,1]]; list.keepY[index[1, 2]];
        #   #cor.upred[index]
        #   summary$u[[k]] = c(list.keepX[index[1,1]], list.keepY[index[1, 2]], cor.upred[index])
        # 
        #   # t.pred
        #   index = which(cor.tpred == max(cor.tpred), arr.ind = TRUE)
        #   #list.keepX[index[1,1]]; list.keepY[index[1, 2]];
        #   #cor.tpred[index]
        #   summary$t[[k]] = c(list.keepX[index[1,1]], list.keepY[index[1, 2]], cor.tpred[index])
        # 
        #   # RSS.u: same as upred
        #   index = which(RSS.upred == min(RSS.upred), arr.ind = TRUE)
        #   #list.keepX[index[1,1]]; list.keepY[index[1, 2]];
        #   #RSS.upred[index]
        #   summary$RSSu[[k]] = c(list.keepX[index[1,1]], list.keepY[index[1, 2]], RSS.upred[index])
        # 
        # 
        #   # RSS.t: different from tpred?
        #   index = which(RSS.tpred == min(RSS.tpred), arr.ind = TRUE)
        #   #list.keepX[index[1,1]]; list.keepY[index[1, 2]];
        #   #RSS.tpred[index]
        #   summary$RSSt[[k]] = c(list.keepX[index[1,1]], list.keepY[index[1, 2]], RSS.tpred[index])
        # } # end ncomp
        #   out$summary = summary
        
        if(isTRUE(pls.model)){
            out$cor.pred = cor.pred
            out$RSS.pred = RSS.pred
            out$Q2.tot.ave = apply(Q2.tot.ave, 1, mean)
        }else{
            if(type.tune == 'cor') out$cor.pred = cor.pred
            if(type.tune == 'RSS') out$RSS.pred = RSS.pred
        }
        
        return(out)
    } 
