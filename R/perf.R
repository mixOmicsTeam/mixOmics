#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#
# created: 2015
# last modified: 05-10-2017
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


# --------------------------
# declare the S3 function:
# -------------------------






#' Compute evaluation criteria for PLS, sPLS, PLS-DA, sPLS-DA, MINT and DIABLO
#'
#' Function to evaluate the performance of the fitted PLS, sparse PLS, PLS-DA,
#' sparse PLS-DA, MINT (mint.splsda) and DIABLO (block.splsda) models using
#' various criteria.
#'
#' Procedure. The process of evaluating the performance of a fitted model
#' \code{object} is similar for all PLS-derived methods; a cross-validation
#' approach is used to fit the method of \code{object} on \code{folds-1}
#' subsets of the data and then to predict on the subset left out. Different
#' measures of performance are available depending on the model. Parameters
#' such as \code{logratio}, \code{multilevel}, \code{keepX} or \code{keepY} are
#' retrieved from \code{object}.
#'
#' Parameters. If \code{validation = "Mfold"}, M-fold cross-validation is
#' performed. \code{folds} specifies the number of folds to generate. The folds
#' also can be supplied as a list of vectors containing the indexes defining
#' each fold as produced by \code{split}. When using \code{validation =
#' "Mfold"}, make sure that you repeat the process several times (as the
#' results will be highly dependent on the random splits and the sample size).
#'
#' If \code{validation = "loo"}, leave-one-out cross-validation is performed
#' (in that case, there is no need to repeat the process).
#'
#' Measures of performance. For fitted PLS and sPLS regression models,
#' \code{perf} estimates the mean squared error of prediction (MSEP),
#' \eqn{R^2}, and \eqn{Q^2} to assess the predictive perfity of the model using
#' M-fold or leave-one-out cross-validation. Note that only the \code{classic},
#' \code{regression} and \code{invariant} modes can be applied. For sPLS, the
#' MSEP, \eqn{R^2}, and \eqn{Q^2} criteria are averaged across all folds. Note
#' that for PLS and sPLS objects, perf is performed on the pre-processed data
#' after log ratio transform and multilevel analysis, if any.
#'
#' Sparse methods. The sPLS, sPLS-DA and sgccda functions are run on several
#' and different subsets of data (the cross-folds) and will certainly lead to
#' different subset of selected features. Those are summarised in the output
#' \code{features$stable} (see output Value below) to assess how often the
#' variables are selected across all folds. Note that for PLS-DA and sPLS-DA
#' objects, perf is performed on the original data, i.e. before the
#' pre-processing step of the log ratio transform and multilevel analysis, if
#' any. In addition for these methods, the classification error rate is
#' averaged across all folds.
#'
#' The mint.sPLS-DA function estimates errors based on Leave-one-group-out
#' cross validation (where each levels of object$study is left out (and
#' predicted) once) and provides study-specific outputs
#' (\code{study.specific.error}) as well as global outputs
#' (\code{global.error}).
#'
#' AUROC. For PLS-DA, sPLS-DA, mint.PLS-DA and mint.sPLS-DA methods: if
#' \code{auc=TRUE}, Area Under the Curve (AUC) values are calculated from the
#' predicted scores obtained from the \code{predict} function applied to the
#' internal test sets in the cross-validation process, either for all samples
#' or for study-specific samples (for mint models). Therefore we minimise the
#' risk of overfitting. See \code{\link{auroc}} for more details. Our
#' multivariate supervised methods already use a prediction threshold based on
#' distances (see \code{predict}) that optimally determine class membership of
#' the samples tested. As such AUC and ROC are not needed to estimate the
#' performance of the model. We provide those outputs as complementary
#' performance measures. See more details in our mixOmics article.
#'
#' Prediction distances. See details from \code{?predict}, and also our
#' supplemental material in the mixOmics article.
#'
#' Repeats of the CV-folds. Repeated cross-validation implies that the whole CV
#' process is repeated a number of times (\code{nrepeat}) to reduce variability
#' across the different subset partitions. In the case of Leave-One-Out CV
#' (\code{validation = 'loo'}), each sample is left out once (\code{folds = N}
#' is set internally) and therefore nrepeat is by default 1.
#'
#' BER is appropriate in case of an unbalanced number of samples per class as
#' it calculates the average proportion of wrongly classified samples in each
#' class, weighted by the number of samples in each class. BER is less biased
#' towards majority classes during the performance assessment.
#'
#' More details about the PLS modes in \code{?pls}.
#'
#' @aliases perf perf.mixo_pls perf.mixo_spls perf.mixo_plsda perf.mixo_splsda
#' perf.mint.splsda perf.sgccda
#' @param object object of class inherited from \code{"pls"}, \code{"plsda"},
#' \code{"spls"}, \code{"splsda"} or \code{"mint.splsda"}. The function will
#' retrieve some key parameters stored in that object.
#' @param dist only applies to an object inheriting from \code{"plsda"},
#' \code{"splsda"} or \code{"mint.splsda"} to evaluate the classification
#' performance of the model. Should be a subset of \code{"max.dist"},
#' \code{"centroids.dist"}, \code{"mahalanobis.dist"}. Default is \code{"all"}.
#' See \code{\link{predict}}.
#' @param validation character.  What kind of (internal) validation to use,
#' matching one of \code{"Mfold"} or \code{"loo"} (see below). Default is
#' \code{"Mfold"}.
#' @param folds the folds in the Mfold cross-validation. See Details.
#' @param nrepeat Number of times the Cross-Validation process is repeated.
#' This is an important argument to ensure the estimation of the performance to
#' be as accurate as possible.
#' @param auc if \code{TRUE} calculate the Area Under the Curve (AUC)
#' performance of the model.
#' @param progressBar by default set to \code{TRUE} to output the progress bar
#' of the computation.
#' @param cpus Number of cpus to use when running the code in parallel.
#' @param ... not used
#' @return For PLS and sPLS models, \code{perf} produces a list with the
#' following components: \item{MSEP}{Mean Square Error Prediction for each
#' \eqn{Y} variable, only applies to object inherited from \code{"pls"}, and
#' \code{"spls"}.} \item{R2}{a matrix of \eqn{R^2} values of the
#' \eqn{Y}-variables for models with \eqn{1, \ldots ,}\code{ncomp} components,
#' only applies to object inherited from \code{"pls"}, and \code{"spls"}.}
#' \item{Q2}{if \eqn{Y} containts one variable, a vector of \eqn{Q^2} values
#' else a list with a matrix of \eqn{Q^2} values for each \eqn{Y}-variable.
#' Note that in the specific case of an sPLS model, it is better to have a look
#' at the Q2.total criterion, only applies to object inherited from
#' \code{"pls"}, and \code{"spls"}} \item{Q2.total}{a vector of \eqn{Q^2}-total
#' values for models with \eqn{1, \ldots ,}\code{ncomp} components, only
#' applies to object inherited from \code{"pls"}, and \code{"spls"}}
#' \item{features}{a list of features selected across the folds
#' (\code{$stable.X} and \code{$stable.Y}) for the \code{keepX} and
#' \code{keepY} parameters from the input object.} \item{error.rate}{ For
#' PLS-DA and sPLS-DA models, \code{perf} produces a matrix of classification
#' error rate estimation. The dimensions correspond to the components in the
#' model and to the prediction method used, respectively. Note that error rates
#' reported in any component include the performance of the model in earlier
#' components for the specified \code{keepX} parameters (e.g. error rate
#' reported for component 3 for \code{keepX = 20} already includes the fitted
#' model on components 1 and 2 for \code{keepX = 20}). For more advanced usage
#' of the \code{perf} function, see \url{www.mixomics.org/methods/spls-da/} and
#' consider using the \code{predict} function.} \item{auc}{Averaged AUC values
#' over the \code{nrepeat}}
#'
#' For mint.splsda models, \code{perf} produces the following outputs:
#' \item{study.specific.error}{A list that gives BER, overall error rate and
#' error rate per class, for each study} \item{global.error}{A list that gives
#' BER, overall error rate and error rate per class for all samples}
#' \item{predict}{A list of length \code{ncomp} that produces the predicted
#' values of each sample for each class} \item{class}{A list which gives the
#' predicted class of each sample for each \code{dist} and each of the
#' \code{ncomp} components. Directly obtained from the \code{predict} output.}
#' \item{auc}{AUC values} \item{auc.study}{AUC values for each study}
#'
#' For sgccda models, \code{perf} produces the following outputs:
#' \item{error.rate}{Prediction error rate for each block of \code{object$X}
#' and each \code{dist}} \item{error.rate.per.class}{Prediction error rate for
#' each block of \code{object$X}, each \code{dist} and each class}
#' \item{predict}{Predicted values of each sample for each class, each block
#' and each component} \item{class}{Predicted class of each sample for each
#' block, each \code{dist}, each component and each nrepeat} \item{features}{a
#' list of features selected across the folds (\code{$stable.X} and
#' \code{$stable.Y}) for the \code{keepX} and \code{keepY} parameters from the
#' input object.} \item{AveragedPredict.class}{if more than one block, returns
#' the average predicted class over the blocks (averaged of the \code{Predict}
#' output and prediction using the \code{max.dist} distance)}
#' \item{AveragedPredict.error.rate}{if more than one block, returns the
#' average predicted error rate over the blocks (using the
#' \code{AveragedPredict.class} output)} \item{WeightedPredict.class}{if more
#' than one block, returns the weighted predicted class over the blocks
#' (weighted average of the \code{Predict} output and prediction using the
#' \code{max.dist} distance)} \item{WeightedPredict.error.rate}{if more than
#' one block, returns the weighted average predicted error rate over the blocks
#' (using the \code{WeightedPredict.class} output)} \item{MajorityVote}{if more
#' than one block, returns the majority class over the blocks. NA for a sample
#' means that there is no consensus on the predicted class for this particular
#' sample over the blocks.} \item{MajorityVote.error.rate}{if more than one
#' block, returns the error rate of the \code{MajorityVote} output}
#' \item{WeightedVote}{if more than one block, returns the weighted majority
#' class over the blocks. NA for a sample means that there is no consensus on
#' the predicted class for this particular sample over the blocks.}
#' \item{WeightedVote.error.rate}{if more than one block, returns the error
#' rate of the \code{WeightedVote} output} \item{weights}{Returns the weights
#' of each block used for the weighted predictions, for each nrepeat and each
#' fold} \item{choice.ncomp}{For supervised models; returns the optimal number
#' of components for the model for each prediction distance using one-sided
#' t-tests that test for a significant difference in the mean error rate (gain
#' in prediction) when components are added to the model. See more details in
#' Rohart et al 2017 Suppl. For more than one block, an optimal ncomp is
#' returned for each prediction framework.}
#' @author Ignacio González, Amrit Singh, Kim-Anh Lê Cao, Benoit Gautier,
#' Florian Rohart.
#' @seealso \code{\link{predict}}, \code{\link{nipals}},
#' \code{\link{plot.perf}}, \code{\link{auroc}} and \url{www.mixOmics.org} for
#' more details.
#' @references DIABLO:
#'
#' Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Lê
#' Cao K.A. (2016). DIABLO - multi omics integration for biomarker discovery.
#'
#' mixOmics article:
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#'
#' MINT:
#'
#' Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017). MINT: A
#' multivariate integrative approach to identify a reproducible biomarker
#' signature across multiple experiments and platforms. BMC Bioinformatics
#' 18:128.
#'
#' PLS and PLS citeria for PLS regression: Tenenhaus, M. (1998). \emph{La
#' regression PLS: theorie et pratique}. Paris: Editions Technic.
#'
#' Chavent, Marie and Patouille, Brigitte (2003). Calcul des coefficients de
#' regression et du PRESS en regression PLS1. \emph{Modulad n}, \bold{30} 1-11.
#' (this is the formula we use to calculate the Q2 in perf.pls and perf.spls)
#'
#' Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction
#' (MSEP) Estimates for Principal Component Regression (PCR) and Partial Least
#' Squares Regression (PLSR). \emph{Journal of Chemometrics} \bold{18}(9),
#' 422-429.
#'
#' sparse PLS regression mode:
#'
#' Lê Cao, K. A., Rossouw D., Robert-Granie, C. and Besse, P. (2008). A sparse
#' PLS for variable selection when integrating Omics data. \emph{Statistical
#' Applications in Genetics and Molecular Biology} \bold{7}, article 35.
#'
#' One-sided t-tests (suppl material):
#'
#' Rohart F, Mason EA, Matigian N, Mosbergen R, Korn O, Chen T, Butcher S,
#' Patel J, Atkinson K, Khosrotehrani K, Fisk NM, Lê Cao K-A&, Wells CA&
#' (2016). A Molecular Classification of Human Mesenchymal Stromal Cells. PeerJ
#' 4:e1845.
#' @keywords regression multivariate
#' @examples
#'
#'
#' ## validation for objects of class 'pls' (regression)
#' # ----------------------------------------
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#'
#'
#' # try tune the number of component to choose
#' # ---------------------
#' # first learn the full model
#' liver.pls <- pls(X, Y, ncomp = 10)
#'
#' # with 5-fold cross validation: we use the same parameters as in model above
#' # but we perform cross validation to compute the MSEP, Q2 and R2 criteria
#' # ---------------------------
#' liver.val <- perf(liver.pls, validation = "Mfold", folds = 5)
#'
#' # Q2 total should decrease until it reaches a threshold
#' liver.val$Q2.total
#'
#' # ncomp = 2 is enough
#' plot(liver.val$Q2.total, type = 'l', col = 'red', ylim = c(-0.5, 0.5),
#' xlab = 'PLS components', ylab = 'Q2 total')
#' abline(h = 0.0975, col = 'darkgreen')
#' legend('topright', col = c('red', 'darkgreen'),
#' legend = c('Q2 total', 'threshold 0.0975'), lty = 1)
#' title('Liver toxicity PLS 5-fold, Q2 total values')
#'
#'
#' \dontrun{
#' #have a look at the other criteria
#' # ----------------------
#' # R2
#' liver.val$R2
#' matplot(t(liver.val$R2), type = 'l', xlab = 'PLS components', ylab = 'R2 for each variable')
#' title('Liver toxicity PLS 5-fold, R2 values')
#'
#' # MSEP
#' liver.val$MSEP
#' matplot(t(liver.val$MSEP), type = 'l', xlab = 'PLS components', ylab = 'MSEP for each variable')
#' title('Liver toxicity PLS 5-fold, MSEP values')
#'
#' ## validation for objects of class 'spls' (regression)
#' # ----------------------------------------
#' ncomp = 7
#' # first, learn the model on the whole data set
#' model.spls = spls(X, Y, ncomp = ncomp, mode = 'regression',
#' keepX = c(rep(10, ncomp)), keepY = c(rep(4,ncomp)))
#'
#'
#' # with leave-one-out cross validation
#' ##set.seed(45)
#' model.spls.val <- perf(model.spls, validation = "Mfold", folds = 5 )#validation = "loo")
#'
#' #Q2 total
#' model.spls.val$Q2.total
#'
#' # R2:we can see how the performance degrades when ncomp increases
#' model.spls.val$R2
#' plot(model.spls.val, criterion="R2", type = 'l')
#' plot(model.spls.val, criterion="Q2", type = 'l')
#'
#'
#' ## validation for objects of class 'splsda' (classification)
#' # ----------------------------------------
#' X <- srbct$gene
#' Y <- srbct$class
#'
#' ncomp = 2
#'
#' srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))
#'
#' # with Mfold
#' # ---------
#' set.seed(45)
#' error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
#' dist = "all", auc = TRUE)
#' error
#' error$auc
#'
#' plot(error)
#'
#' # parallel code
#' set.seed(45)
#' error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
#' dist = "all", auc = TRUE, cpus =2)
#'
#'
#'
#'
#'
#' # with 5 components and nrepeat =5, to get a $choice.ncomp
#' ncomp = 5
#' srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))
#'
#' set.seed(45)
#' error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
#' dist = "all", nrepeat =5)
#' error
#'
#' plot(error)
#'
#' # parallel code
#' set.seed(45)
#' error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
#' dist = "all", auc = TRUE, cpus =2)
#'
#'
#'
#' ## validation for objects of class 'mint.splsda' (classification)
#' # ----------------------------------------
#'
#' res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3, keepX = c(10, 5, 15),
#' study = stemcells$study)
#'
#' out = perf(res, auc = TRUE)
#' out
#'
#' out$auc
#' out$auc.study
#'
#' ## validation for objects of class 'sgccda' (classification)
#' # ----------------------------------------
#'
#' Y = nutrimouse$diet
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
#' design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#'
#' nutrimouse.sgccda <- block.splsda(X=data,
#' Y = Y,
#' design = design,
#' keepX = list(gene=c(10,10), lipid=c(15,15)),
#' ncomp = 2,
#' scheme = "horst")
#'
#' perf = perf(nutrimouse.sgccda)
#' perf
#'
#'
#'
#' #with 5 components and nrepeat=5 to get $choice.ncomp
#' nutrimouse.sgccda <- block.splsda(X=data,
#' Y = Y,
#' design = design,
#' keepX = list(gene=c(10,10), lipid=c(15,15)),
#' ncomp = 5,
#' scheme = "horst")
#'
#' perf = perf(nutrimouse.sgccda, folds = 5, nrepeat = 5)
#' perf
#'
#' perf$choice.ncomp
#'
#' }
#'
#' @export perf
perf = function(object, ...) UseMethod("perf")


#------------------------------------------------------#
#-- Includes perf for PLS, sPLS, PLS-DA and sPLS-DA --#
#------------------------------------------------------#

#---------------------------------------------------
# perf for spls and pls object
#---------------------------------------------------

perf.mixo_spls  = perf.mixo_pls = function(object,
validation = c("Mfold", "loo"),
folds = 10,
progressBar = TRUE,
...)
{
    #------------------#
    #-- check entries --#

    #-- check spls mode
    if (object$mode == 'canonical')
    stop("object$mode should be 'regression', 'invariant' or 'classic'.", call. = FALSE)

    #-- validation
    choices = c("Mfold", "loo")
    validation = choices[pmatch(validation, choices)]

    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)

    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)

    #-- end checking --#
    #------------------#


    #-- cross-validation approach  ---------------------------------------------#
    #---------------------------------------------------------------------------#

    #-- initialising arguments --#
    # these are the centered and scaled matrices output from pls, we remove $nzv if needed
    if (length(object$nzv$Position)>0)
    {
        X = object$X[, -object$nzv$Position]
    } else {
        X = object$X
    }
    Y = object$Y

    scale = object$scale
    tol = object$tol
    max.iter = object$max.iter
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()

    if (any(is.na(X)) || any(is.na(Y)))
    stop("missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.", call. = FALSE)


    #-- tells which variables are selected in X and in Y --#
    if (is(object, "mixo_spls"))
    {
        keepX = object$keepX
        keepY = object$keepY
    } else {
        keepX = rep(ncol(X), ncomp)
        keepY = rep(ncol(Y), ncomp)
    }

    #-- define the folds --#
    if (validation == "Mfold")
    {
        if (is.list(folds))
        {

            if (length(folds) < 2 || length(folds) > n)
            stop("Invalid number of folds.", call. = FALSE)

            if (length(unlist(folds)) != n)
            stop("Invalid folds. The total number of samples in folds must be equal to ",
            n, ".", call. = FALSE)

            if (length(unique(unlist(folds))) != n)
            stop("Invalid folds. Repeated samples in folds.", call. = FALSE)

            M = length(folds)
        } else {
            if (is.null(folds) || !is.finite(folds) || folds < 2 || folds > n)
            {
                stop("Invalid number of folds.", call. = FALSE)
            } else {
                M = round(folds)
                folds = split(sample(1:n), rep(1:M, length = n))
            }
        }
    } else {
        folds = split(1:n, rep(1:n, length = n))
        M = n
    }

    #-- set up a progress bar --#
    if (progressBar == TRUE)
    {
        pb = txtProgressBar(style = 3)
        nBar = 1
    }

    #-- initialize new objects --#
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    PRESS.inside = Q2 = MSEP = R2 = matrix(nrow = ncomp, ncol = q)
    MSEP.mat = Ypred = array(0, c(n, q, ncomp))

    press.mat = lapply(1 : ncomp, function(x){matrix(NA, nrow = n, ncol = q)})
    RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = q)})
    RSS.indiv[[1]] = X

    #-- record feature stability --#
    # initialize new objects:= to record feature stability
    featuresX  = featuresY =  list()
    for(k in 1:ncomp)
    featuresX[[k]] = featuresY[[k]] = NA


    #-- loop on h = ncomp --#
    # the loop is only for the calculation of Q2 on each component
    for (h in 1:ncomp)
    {

        #-- initialising arguments --#
        tt = object$variates$X[, h]
        u = object$variates$Y[, h]
        b = object$loadings$Y[, h]
        nx = p - keepX[h]
        ny = q - keepY[h]

        #only used for matrices deflation
        c = crossprod(X, tt)/drop(crossprod(tt)) #object$mat.c[, h]
        d = crossprod(Y, tt)/drop(crossprod(tt))#object$mat.d[, h]

        RSS.indiv[[h + 1]] = Y - tt %*% t(d)
        RSS[h + 1, ] = colSums((Y - tt %*% t(d))^2)

        #-- loop on i (cross validation) --#
        for (i in 1:M)
        {
            if (progressBar == TRUE)
            {
                setTxtProgressBar(pb, nBar/(ncomp * M))
                nBar = nBar + 1
            }

            omit = folds[[i]]
            X.train = X[-omit, , drop = FALSE]
            Y.train = Y[-omit, , drop = FALSE]
            X.test = X[omit, , drop = FALSE]
            Y.test = Y[omit, , drop = FALSE]
            u.cv = u[-omit]


            #-- for MSEP and R2 criteria, no loop on the component as we do a spls with ncomp
            if (h == 1)
            {
                #nzv = (apply(X.train, 2, var) > .Machine$double.eps) # removed in v6.0.0 so that MSEP, R2 and Q2 are obtained with the same data
                # re-added in >6.1.3 to remove constant variables
                nzv = (apply(X.train, 2, var) > .Machine$double.eps)

                # creating a keepX.temp that can change for each fold, depending on nzv
                keepX.temp = keepX
                if(any(keepX.temp > sum(nzv)))
                keepX.temp[which(keepX.temp>sum(nzv))] = sum(nzv)

                spls.res = mixOmics::spls(X.train[,nzv], Y.train, ncomp = ncomp, mode = mode, max.iter = max.iter, tol = tol, keepX = keepX.temp, keepY = keepY, near.zero.var = FALSE, scale = scale)
                Y.hat = predict.mixo_spls(spls.res, X.test[,nzv, drop = FALSE])$predict
                if(sum(is.na(Y.hat))>0) break
                for (k in 1:ncomp)
                {
                    Ypred[omit, , k] = Y.hat[, , k]
                    MSEP.mat[omit, , k] = (Y.test - Y.hat[, , k])^2

                    # added: record selected features in each set
                    if (is(object,"mixo_spls"))
                    {
                        featuresX[[k]] = c(unlist(featuresX[[k]]), selectVar(spls.res, comp = k)$X$name)
                        featuresY[[k]] = c(unlist(featuresY[[k]]), selectVar(spls.res, comp = k)$Y$name)
                    }

                } # end loop on k
            }

            #-- Q2 criterion
            a.old.cv = 0
            iter.cv = 1

            repeat{
                a.cv = crossprod(X.train, u.cv)
                if (nx != 0)
                {
                    a.cv = ifelse(abs(a.cv) > abs(a.cv[order(abs(a.cv))][nx]),
                    (abs(a.cv) - abs(a.cv[order(abs(a.cv))][nx])) * sign(a.cv), 0)
                }
                a.cv = a.cv / drop(sqrt(crossprod(a.cv)))
                t.cv = X.train %*% a.cv

                b.cv = crossprod(Y.train, t.cv)
                if (ny != 0)
                {
                    b.cv = ifelse(abs(b.cv) > abs(b.cv[order(abs(b.cv))][ny]),
                    (abs(b.cv) - abs(b.cv[order(abs(b.cv))][ny])) * sign(b.cv), 0)
                }
                b.cv = b.cv / drop(sqrt(crossprod(b.cv)))
                u.cv = Y.train %*% b.cv

                if ((crossprod(a.cv - a.old.cv) < tol) || (iter.cv == max.iter))
                break

                a.old.cv = a.cv
                iter.cv = iter.cv + 1
            }

            d.cv = t(Y.train) %*% (X.train %*% a.cv) / norm((X.train %*% a.cv), type = "2")^2
            Y.hat.cv = (X.test %*% a.cv) %*% t(d.cv)
            press.mat[[h]][omit, ] = Y.test - Y.hat.cv


        } # end i (cross validation)

        #-- compute the Q2 creterion --#
        PRESS.inside[h, ] = apply(press.mat[[h]], 2, function(x){norm(x, type = "2")^2})
        Q2[h, ] = 1 - PRESS.inside[h, ] / RSS[h, ]

        #-- deflation des matrices (for Q2 criterion)
        X = X - tt %*% t(c)

        #-- mode classic
        if (mode == "classic")
        Y = Y - tt %*% t(b)

        #-- mode regression
        if (mode == "regression")
        Y = Y - tt %*% t(d)

        #-- mode invariant: Y is unchanged

        #-- compute the MSEP creterion --#
        MSEP[h, ] = apply(as.matrix(MSEP.mat[, , h]), 2, mean)

        #-- compute the R2 creterion --#
        R2[h, ] = (diag(cor(object$Y, Ypred[, , h])))^2

    } #-- end loop on h --#

    if (progressBar == TRUE) cat('\n')


    #-- output -----------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    Q2.total = matrix(1 - rowSums(PRESS.inside) / rowSums(RSS[-(ncomp+1), , drop = FALSE]), nrow = 1, ncol = ncomp,
    dimnames = list("Q2.total", paste0(1:ncomp, " comp")))

    # set up dimnames
    rownames(MSEP) = rownames(R2) = rownames(Q2) = paste0(1:ncomp, " comp")
    colnames(MSEP) = colnames(R2) = colnames(Q2) = object$names$colnames$Y

    result = list()
    result$MSEP = t(MSEP)
    result$R2 = t(R2)
    result$Q2 = t(Q2)
   	result$Q2.total =  t(Q2.total)
    result$RSS = RSS
    result$PRESS = PRESS.inside
    result$press.mat = press.mat
    result$RSS.indiv = RSS.indiv

    #---- extract stability of features -----#
    if (is(object, "mixo_spls"))
    {
        list.features.X = list()
        list.features.Y = list()

        for(k in 1:ncomp)
        {
            #remove the NA value that was added for initialisation
            remove.naX = which(is.na(featuresX[[k]]))
            remove.naY = which(is.na(featuresY[[k]]))
            # then summarise as a factor and output the percentage of appearance
            list.features.X[[k]] = sort(table(as.factor(featuresX[[k]][-remove.naX])) / M, decreasing = TRUE)
            list.features.Y[[k]] = sort(table(as.factor(featuresY[[k]][-remove.naY])) / M, decreasing = TRUE)

        }
        names(list.features.X)  = names(list.features.Y) = paste0('comp', 1:ncomp)

        # features
        result$features$stable.X = list.features.X
        result$features$stable.Y = list.features.Y
    }

    #--- class
    if (is(object,"mixo_spls"))
    {
        method = "spls.mthd"
    } else if (is(object, "mixo_pls")) {
        method = "pls.mthd"
    } else {
        warning("Something that should not happen happened. Please contact us.")
    }
    class(result) = c("perf",paste(c("perf", method), collapse ="."))
    result$call = match.call()

    return(invisible(result))
}


# ---------------------------------------------------
# perf for plsda and splsda object
# ---------------------------------------------------
perf.mixo_splsda = perf.mixo_plsda = function(object,
dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
validation = c("Mfold", "loo"),
folds = 10,
nrepeat = 1,
auc = FALSE,
progressBar = TRUE,
cpus,
...)
{

    #-- initialising arguments --#
    # these data are the centered and scaled X output or the unmapped(Y) scaled and centered
    X = object$input.X
    level.Y = object$names$colnames$Y  #to make sure the levels are ordered
    Y = object$Y
    ncomp = object$ncomp
    n = nrow(X)

    logratio = object$logratio
    if (is.null(logratio))
    logratio = "none"

    multilevel = object$multilevel # repeated measurement and Y
    near.zero.var = !is.null(object$nzv) # if near.zero.var was used, we set it to TRUE. if not used, object$nzv is NULL

    #-- tells which variables are selected in X and in Y --#

    if (is(object, "mixo_splsda"))
    {
        keepX = object$keepX
    } else {
        keepX = rep(ncol(X), ncomp)
    }

    tol = object$tol
    max.iter = object$max.iter
    scale = object$scale

    # initialize new objects:
    features = list()
    for(k in 1:ncomp)
    features[[k]] = NA

    # check input arguments

    if (hasArg(method.predict))
    stop("'method.predict' argument has been replaced by 'dist' to match the 'tune' function")
    method.predict = NULL # to pass R CMD check

    dist = match.arg(dist, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
    if (any(dist == "all"))
    {
        nmthdd = 3
        dist = c("max.dist", "centroids.dist", "mahalanobis.dist")
    } else {
        nmthdd = length(dist)
    }

    if (length(validation) > 1 )
    validation = validation [1]
    if (!(validation %in% c("Mfold", "loo")))
    stop("Choose 'validation' among the two following possibilities: 'Mfold' or 'loo'")

    if (validation == "loo")
    {
        if (nrepeat != 1)
        warning("Leave-One-Out validation does not need to be repeated: 'nrepeat' is set to '1'.")
        nrepeat = 1
    }

    if (!is.logical(progressBar))
    stop("'progressBar' must be either TRUE or FALSE")

    measure = c("overall","BER") # one of c("overall","BER")


    if (!(logratio %in% c("none", "CLR")))
    stop("Choose one of the two following logratio transformation: 'none' or 'CLR'")
    #fold is checked in 'MCVfold'

    if (!isNULL(cpus))
    {
        if(!is.numeric(cpus) | length(cpus)!=1)
        stop("'cpus' must be a numerical value")

        parallel = TRUE
        cl = makeCluster(cpus, type = "SOCK")
        #clusterExport(cl, c("splsda","selectVar"))
        clusterEvalQ(cl, library(mixOmics))

    } else {
        parallel = FALSE
        cl = NULL
    }


    #---------------------------------------------------------------------------#
    #-- logration + multilevel approach ----------------------------------------#
    # we can do logratio and multilevel on the whole data as these transformation are done per sample
    X = logratio.transfo(X = X, logratio = logratio)
    if (!is.null(multilevel))
    {
        Xw = withinVariation(X, design = multilevel)
        X = Xw
    }
    #-- logratio + multilevel approach -----------------------------------------#
    #---------------------------------------------------------------------------#


    # -------------------------------------
    # added: first check for near zero var on the whole data set
    if (near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
            X = X[, -nzv$Position, drop=TRUE]

            if (ncol(X)==0)
            stop("No more predictors after Near Zero Var has been applied!")

            if (any(keepX > ncol(X)))
            keepX = ncol(X)

        }
    }
    # and then we start from the X data set with the nzv removed

    #---------------------------------------------------------------------------#
    #-- NA calculation      ----------------------------------------------------#

    misdata = c(X=anyNA(X), Y=FALSE) # Detection of missing data. we assume no missing values in the factor Y

    if (any(misdata))
    {
        is.na.A = is.na(X)

        #ind.NA = which(apply(is.na.A, 1, sum) > 0) # calculated only once
        #ind.NA.col = which(apply(is.na.A, 2, sum) > 0) # calculated only once
    } else {
        is.na.A = NULL
        #ind.NA = ind.NA.col = NULL
    }
    #-- NA calculation      ----------------------------------------------------#
    #---------------------------------------------------------------------------#


    list.features = list()

    mat.error.rate = mat.sd.error = mat.mean.error = error.per.class.keepX.opt = error.per.class.keepX.opt.mean = list()
    error.per.class = list()
    final=list()

    for (measure_i in measure)
    {
        mat.sd.error[[measure_i]] = matrix(0,nrow = ncomp, ncol = length(dist),
        dimnames = list(c(paste0('comp', 1 : ncomp)), dist))
        mat.mean.error[[measure_i]] = matrix(0,nrow = ncomp, ncol = length(dist),
        dimnames = list(c(paste0('comp', 1 : ncomp)), dist))
        error.per.class.keepX.opt[[measure_i]] = list()
        error.per.class.keepX.opt.mean[[measure_i]] = list()
        mat.error.rate[[measure_i]]=list()
        for(ijk in dist)
        {
            mat.error.rate[[measure_i]][[ijk]] = matrix(0, nrow = ncomp, ncol = nrepeat,
            dimnames = list(c(paste0('comp', 1 : ncomp)), c(paste0('nrep', 1 : nrepeat))))

            error.per.class.keepX.opt[[measure_i]][[ijk]] = array(0, c(nlevels(Y), nrepeat, ncomp),
            dimnames = list(c(levels(Y)), c(paste0('nrep', 1 : nrepeat)), c(paste0('comp', 1:ncomp, sep=''))))

            error.per.class.keepX.opt.mean[[measure_i]][[ijk]] = matrix(nrow = nlevels(Y), ncol = ncomp,
            dimnames = list(c(levels(Y)), c(paste0('comp', 1 : ncomp))))
        }
    }

    if(auc == TRUE)
    {
        auc.mean=list()
        auc.all=list()
    }

    prediction.all = class.all = auc.mean = auc.all = list()
    for(ijk in dist)
    {
        class.all[[ijk]] = array(0, c(nrow(X),  nrepeat ,ncomp),
        dimnames = list(rownames(X),c(paste0('nrep', 1 : nrepeat)),c(paste0('comp', 1 : ncomp))))
    }

    class.object=class(object)
    if (!isNULL(cpus))
    clusterExport(cl, c("X","Y","is.na.A","misdata","scale","near.zero.var","class.object","test.keepX"),envir=environment())

    for (comp in 1 : ncomp)
    {
        if (progressBar == TRUE)
        cat("\ncomp",comp, "\n")


        if(comp > 1)
        {
            choice.keepX = keepX[1 : (comp - 1)]
        } else {
           choice.keepX = NULL
        }
        test.keepX = keepX[comp]
        names(test.keepX) = test.keepX
        #test.keepX is a value

        # estimate performance of the model for each component
        result = MCVfold.spls (X, Y, multilevel = multilevel, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = comp,
        choice.keepX = choice.keepX, test.keepX = test.keepX, test.keepY = nlevels(Y),
        measure = measure, dist = dist, scale=scale,
        near.zero.var = near.zero.var,
        auc = auc, progressBar = progressBar, class.object = class.object, cl = cl, parallel = parallel,
        misdata = misdata, is.na.A = is.na.A)#, ind.NA = ind.NA, ind.NA.col = ind.NA.col)

        # ---- extract stability of features ----- # NEW
        if (is(object, "mixo_splsda"))
        list.features[[comp]] = result$features$stable

        for (ijk in dist)
        {
            for (measure_i in measure)
            {
                mat.error.rate[[measure_i]][[ijk]][comp,] = result[[measure_i]]$mat.error.rate[[ijk]][1,]
                mat.mean.error[[measure_i]][comp, ijk]=result[[measure_i]]$error.rate.mean[[ijk]]
                if (!is.null(result[[measure_i]]$error.rate.sd))
                {
                    mat.sd.error[[measure_i]][comp, ijk]=result[[measure_i]]$error.rate.sd[[ijk]]
                } else {
                    mat.sd.error= NULL
                }
                # confusion matrix for keepX.opt, for each nrep
                error.per.class.keepX.opt[[measure_i]][[ijk]][ , ,comp] = result[[measure_i]]$confusion[[ijk]]

                # confusion matrix for keepX.opt, averaged over all nrep
                error.per.class.keepX.opt.mean[[measure_i]][[ijk]][ ,comp] = apply(result[[measure_i]]$confusion[[ijk]],1 , mean)
            }

            #prediction of each samples for each fold and each repeat, on each comp
            class.all[[ijk]][, , comp] = result$class.comp[[ijk]][,,1]
        }
        prediction.all[[comp]] = array(unlist(result$prediction.comp),c(nrow(result$prediction.comp[[1]]), ncol(result$prediction.comp[[1]]), nrepeat),
        dimnames = c(dimnames(result$prediction.comp[[1]])[1:2], list(paste0("nrep",1:nrepeat))))#[[1]][, , 1] #take only one component [[1]] and one of test.keepX [,,1]

        if(auc == TRUE)
        {
            auc.all[[comp]] = lapply(result$auc.all, function(x) x[,,1])
            auc.mean[[comp]] = result$auc[, , 1]
        }
    }
    if (parallel == TRUE)
    stopCluster(cl)

    names(prediction.all) = paste0('comp', 1:ncomp)

    # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
    ncomp_opt = matrix(NA, nrow = length(measure), ncol = length(dist),
    dimnames = list(measure, dist))
    if(nrepeat > 2 & ncomp >1)
    {
        for (measure_i in measure)
        {
            for (ijk in dist)
            ncomp_opt[measure, ijk] = t.test.process(t(mat.error.rate[[measure_i]][[ijk]]))
        }
    }

    result = list(error.rate = mat.mean.error,
    error.rate.sd = mat.sd.error,
    error.rate.all = mat.error.rate,
    error.rate.class = error.per.class.keepX.opt.mean[[1]],
    error.rate.class.all = error.per.class.keepX.opt[[1]],
    predict = prediction.all,
    class = class.all,
    choice.ncomp = ncomp_opt)

    if(auc)
    {
        names(auc.mean) = c(paste0('comp', 1:ncomp))
        result$auc = auc.mean

        names(auc.all) = c(paste0('comp', 1:ncomp))
        result$auc.all =auc.all
    }

    if (is(object, "mixo_splsda"))
    {
        names(list.features) = paste0('comp', 1:ncomp)
        result$features$stable = list.features
    }

    if (progressBar == TRUE)
    cat('\n')

    # added
    if (near.zero.var == TRUE)
    result$nzvX = nzv$Position

    if (is(object, "mixo_splsda"))
    {
        method = "splsda.mthd"
    } else if (is(object, "mixo_plsda")) {
        method = "plsda.mthd"
    } else {
        warning("Something that should not happen happened. Please contact us.")
    }
    class(result) = c("perf",paste(c("perf", method), collapse ="."))
    result$call = match.call()


    #updated outputs
    return(invisible(result))
}


