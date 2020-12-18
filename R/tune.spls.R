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
# tune.spls: chose the optimal number of parameters per component on a spls method
# ========================================================================================================
# TODO details on nrepeat
# TODO details on pls vs spls
# TODO tidy outputs preferrably in an array
#' Tuning functions for sPLS and PLS functions
#' 
#' @template description/tune
#' 
#' @template section/folds
#' @template section/nrepeat
#' @template section/measure-pls
#' @template section/t-test-process
#' 
#' @section more:
#' See also \code{?perf} for more details.
#' 
#' @inheritParams spls
#' @template arg/test.keepX-X.matrix
#' @template arg/test.keepY
#' @template arg/validation
#' @template arg/folds
#' @template arg/nrepeat
#' @param measure One of c('cor', 'RSS') indicating the tuning measure. See
#'   details.
#' @template arg/progressBar
#' @template arg/BPPARAM
#' @param LimQ2 Q2 threshold for recommending optimal \code{ncomp}.
#' @param ... Optional parameters passed to \code{\link{spls}}
#' @return A list that contains: \item{cor.pred}{The correlation of predicted vs
#'   actual components from X (t) and Y (u) for each
#'   component}\item{RSS.pred}{The Residual Sum of Squares of predicted vs
#'   actual components from X (t) and Y (u) for each component}
#'   \item{choice.keepX}{returns the number of variables selected for X (optimal
#'   keepX) on each component.} \item{choice.keepY}{returns the number of
#'   variables selected for Y (optimal keepY) on each component.}
#'   \item{choice.ncomp}{returns the optimal number of components for the model
#'   fitted with \code{$choice.keepX} and \code{$choice.keepY} } \item{call}{The
#'   functioncal call including the parameteres used.}
#' @author Kim-Anh Lê Cao, Al J Abadi, Benoit Gautier, Francois Bartolo,
#' Florian Rohart,
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
#' \dontrun{
#' data(liver.toxicity)
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#' set.seed(42)
#' tune.res = tune.spls( X, Y, ncomp = 3,
#'                   test.keepX = c(5, 10, 15),
#'                   test.keepY = c(3, 6, 8), measure = "cor",
#'                   folds = 5, nrepeat = 3, progressBar = TRUE)
#' tune.res$choice.ncomp
#' tune.res$choice.keepX
#' tune.res$choice.keepY
#' # plot the results
#' plot(tune.res)
#' }
# change this so that it simply wraps perf
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
             measure = c('cor', 'RSS'), ## only if spls model
             BPPARAM = SerialParam(),
             progressBar = FALSE,
             limQ2 = 0.0975,
             ...
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
        
        test.keepX <- .change_if_null(arg = test.keepX, default = ncol(X))
        test.keepY <- .change_if_null(arg = test.keepY, default = ncol(Y))
        
        test.keepX <- unique(test.keepX)
        test.keepY <- unique(test.keepY)
        
        spls.model <- !(length(test.keepX) == 1 & length(test.keepY) == 1)
        
        if (ncol(Y) == 1) {
          res <- tune.spls1(X = X, 
                            Y = Y,
                            ncomp = ncomp,
                            test.keepX = test.keepX,
                            validation = validation,
                            folds = folds,
                            measure = measure, # can do a R2 per Y (correlation, linear regression R2), need to call MSEP (see perf.spls).
                            progressBar = progressBar,
                            nrepeat = nrepeat,
                            ...
          )
            ## --- call
            res$call <- NULL
            ## eval all but X and Y
            mc <- mget(names(formals())[-1:-2], sys.frame(sys.nframe()))
            ## replace function, X and Y with unevaluated call
            mc <- as.call(c(as.list(match.call())[1:3], mc))
            res <- c(list(call = mc), res)
            class(res) <- 'tune.spls1'
            return(res)
        }

        measure <- match.arg(measure, choices = c('cor', 'RSS'))
        
        choice.keepX = choice.keepY = NULL
        measure.table.cols <- c('comp', 'keepX', 'keepY', 'repeat', 't', 'u', 'cor', 'RSS')
        # if (spls.model)
        # {
          measure.pred <- expand.grid(keepX = test.keepX, 
                                      keepY = test.keepY,
                                      V = c('u', 't'),
                                      measure = c('cor', 'RSS'),
                                      comp = seq_len(ncomp),
                                      optimum.keepA = FALSE)
        # } else {
        #   # Q2 only
        #   measure.pred <- expand.grid(
        #                               comp = seq_len(ncomp),
        #                               optimum = FALSE)
        # }
        measure.pred <- data.frame(measure.pred)
        measure.pred <- cbind(measure.pred, 
                              value.v = I(rep(list(data.frame(matrix(NA_real_, ncol= 1, nrow = nrepeat))), 
                                              times = nrow(measure.pred))),
                              value.Q2.total = I(rep(list(data.frame(matrix(NA_real_, ncol= 1, nrow = nrepeat))), 
                                              times = nrow(measure.pred)))
                              )
        # measure.pred$value <- lapply(measure.pred$value, as.data.frame)
        
        use_progressBar <- progressBar & (is(BPPARAM, 'SerialParam'))
        n_keepA <- length(test.keepX) * length(test.keepY)
        
        ## initialise optimal keepX/Y
        measure.pred$optimum.keepA <- FALSE 
        measure.pred[
          measure.pred$keepX == test.keepX[1] &
            measure.pred$keepY == test.keepY[1]
          ,]$optimum.keepA <- TRUE
        
            for (comp in seq_len(ncomp)){
              # TODO tune.pls progressBar should use perf
                if (use_progressBar) {
                    n_tested <- 0
                    cat(sprintf("\ntuning component: %s\n", comp))
                }
                    for(keepX in 1:length(test.keepX)){
                        for(keepY in 1:length(test.keepY)){
                            if (use_progressBar) {
                                n_tested <- n_tested + 1
                                prog_level <- n_tested / n_keepA
                                .progressBar(prog_level, title = 'of features tested')
                            }
                            pls.model <- spls(X = X, Y = Y, 
                                              keepX = c(choice.keepX, test.keepX[keepX]), 
                                              keepY = c(choice.keepY, test.keepY[keepY]), 
                                              ncomp = comp, mode = mode, ...)
                            
                            pls.perf <- perf(pls.model, validation = validation, folds = folds, nrepeat = nrepeat)
                            ## why value.u/t is constant coming out of this?
                            ## now that measure.pred is different for the two, account for it

                                  for (measure_i in c('cor', 'RSS')) ## calculate both but use measure only
                                  {
                                    
                                    for (v in c('u', 't'))
                                    {
                                      ## populate the table for both measures
                                      measure.vpred <- pls.perf$measures[[sprintf("%s.%spred", measure_i, v)]]$values
                                      measure.vpred <- measure.vpred[measure.vpred$comp == comp,]
             
                                      measure.pred[measure.pred$comp == comp & 
                                                     measure.pred$keepX == test.keepX[keepX] &
                                                     measure.pred$keepY == test.keepY[keepY] &
                                                     measure.pred$V == v &
                                                     measure.pred$measure == measure_i
                                                   ,]$value.v<- measure.vpred$value
         
                                      
                                    }
                                    value.Q2.total <- pls.perf$measures$Q2.total$values
                                    value.Q2.total <- filter(value.Q2.total, comp == comp)$value
                                    
                                    measure.pred[measure.pred$comp == comp & 
                                                   measure.pred$keepX == test.keepX[keepX] &
                                                   measure.pred$keepY == test.keepY[keepY] &
                                                   measure.pred$V == 'u' 
                                                 ,]$value.Q2.total <- value.Q2.total
                                    
                                  }

                                  ## optimum only uses measure
                            optimum.u <- measure.pred[measure.pred$comp == comp & 
                                                        measure.pred$optimum.keepA == TRUE &
                                                        measure.pred$V == 'u' &
                                                        measure.pred$measure == measure
                                                      ,]$value.v[[1]]
                            optimum.t <- measure.pred[measure.pred$comp == comp & 
                                                        measure.pred$optimum.keepA == TRUE &
                                                        measure.pred$V == 't' &
                                                        measure.pred$measure == measure
                                                      ,]$value.v[[1]]
                            value.u <- measure.pred[measure.pred$comp == comp & 
                                                      measure.pred$keepX == test.keepX[keepX] &
                                                      measure.pred$keepY == test.keepY[keepY] &
                                                      measure.pred$V == 'u' &
                                                      measure.pred$measure == measure
                                                    ,]$value.v[[1]]
                            value.t <- measure.pred[measure.pred$comp == comp & 
                                                      measure.pred$keepX == test.keepX[keepX] &
                                                      measure.pred$keepY == test.keepY[keepY] &
                                                      measure.pred$V == 't' &
                                                      measure.pred$measure == measure
                                                    ,]$value.v[[1]]
         
                                    ## workaround for constant values in t.test # TODO handle it properly
                                    offset.eps <- seq(1, 2, length.out = length(optimum.u))/1e6
                                    optimum.t <- optimum.t + offset.eps
                                    optimum.u <- optimum.u + offset.eps
                                    value.t <- value.t + offset.eps
                                    value.u <- value.u + offset.eps
                                    
                                    .check_improvement <- function(opt, value, measure, nrepeat) {
                                      ## output TRUE if value improves opt for cor or RSS, else FALSE
                                      if (nrepeat > 2) {
                                        t.test.res <- tryCatch(t.test(x = opt, y = value, alternative = ifelse(measure == 'cor', 'greater', 'less')), error = function(e) e)
                                        improved <- t.test.res$p.value < 0.05
                                        } else
                                        {
                                          ## compare values
                                          improved <-  if (measure == 'cor') mean(opt) < mean(value) else mean(opt) > mean(value)
                                        }
                                      improved
                                    }
                                    
                                    u.improved <-.check_improvement(opt = optimum.u, value = value.u, measure = measure, nrepeat = nrepeat)
                                    t.improved <-.check_improvement(opt = optimum.t, value = value.t, measure = measure, nrepeat = nrepeat)
                                    improved <- if (mode == 'canonical') u.improved & t.improved else u.improved
                                  if (improved)
                                  {
                                      ## set previous to FALSE
                                      measure.pred[measure.pred$comp == comp
                                                   ,]$optimum.keepA <- FALSE
                                      ## update optimum
                                    measure.pred[measure.pred$comp == comp & 
                                                   measure.pred$keepX == test.keepX[keepX] &
                                                   measure.pred$keepY == test.keepY[keepY]
                                                 ,]$optimum.keepA <- TRUE
                                  }
                                
                        } # end keepY
                    } #end keepX
  
              choice.keepX.ncomp <-  measure.pred[measure.pred$comp == comp & 
                                                    measure.pred$optimum.keepA == TRUE &
                                                    measure.pred$measure == measure & 
                                                    measure.pred$V == 't' ## doesn't matter t or u
                                                  ,]$keepX
              choice.keepY.ncomp <-  measure.pred[measure.pred$comp == comp & 
                                                    measure.pred$optimum.keepA == TRUE &
                                                    measure.pred$measure == measure &
                                                    measure.pred$V == 't'
                                                  ,]$keepY
              choice.keepX = c(choice.keepX, choice.keepX.ncomp)
              choice.keepY = c(choice.keepY, choice.keepY.ncomp)
              
            } # end comp
        
        choice.ncomp <- 1
        for (comp in seq_len(ncomp))
        {
          Q2.total <- measure.pred[measure.pred$comp == comp & 
                                     measure.pred$keepX == test.keepX[keepX] &
                                     measure.pred$keepY == test.keepY[keepY] &
                                     measure.pred$V == 'u'
                                   ,]$value.Q2.total[[1]]
          Q2.opt <- measure.pred[measure.pred$comp == choice.ncomp & 
                                     measure.pred$keepX == test.keepX[keepX] &
                                     measure.pred$keepY == test.keepY[keepY] &
                                   measure.pred$V == 'u'
                                   ,]$value.Q2.total[[1]]
          keep.comp <- mean(Q2.total) >= limQ2
          if (keep.comp)
            choice.ncomp <- comp
        }
        names(choice.keepX) <- names(choice.keepY) <- paste0('comp', seq_len(ncomp))
        out$choice.keepX = choice.keepX
        out$choice.keepY = choice.keepY
        out$choice.ncomp = choice.ncomp
        
        ## add mean and sd
        val <- unclass(measure.pred$value.v)
        measure.pred$mean <- sapply(measure.pred$value.v, function(x){
          mean(c(as.matrix(x)), na.rm = TRUE)
        })
        
        measure.pred$sd <- sapply(measure.pred$value.v, function(x){
          sd(c(as.matrix(x)), na.rm = TRUE)
        })
        
        out$measure.pred = measure.pred
        ### evaluate all for output except X and Y to save memory
        ## eval all but X and Y
        mc <- mget(names(formals())[-1:-2], sys.frame(sys.nframe()))
        ## replace function, X and Y with unevaluated call
        mc <- as.call(c(as.list(match.call())[1:3], mc))
        out <- c(list(call = mc), out)
        class <- c('tune.pls')
        class(out) <- if (spls.model) class <- c(class, 'tune.spls')
        return(out)
    } 
